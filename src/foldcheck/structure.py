"""Protein coordinate containers and PDB / mmCIF readers."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from foldcheck.geometry import element_from_name


@dataclass
class Atom:
    name: str
    element: str
    resname: str
    chain: str
    resseq: int
    icode: str
    xyz: np.ndarray
    bfactor: float = 0.0
    occupancy: float = 1.0
    altloc: str = ""

    def copy(self) -> Atom:
        return Atom(
            name=self.name,
            element=self.element,
            resname=self.resname,
            chain=self.chain,
            resseq=self.resseq,
            icode=self.icode,
            xyz=np.array(self.xyz, dtype=float, copy=True),
            bfactor=self.bfactor,
            occupancy=self.occupancy,
            altloc=self.altloc,
        )


@dataclass
class Residue:
    resname: str
    chain: str
    resseq: int
    icode: str
    atoms: dict[str, Atom] = field(default_factory=dict)

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)


@dataclass
class Structure:
    atoms: list[Atom]
    name: str = ""

    def copy(self) -> Structure:
        return Structure(atoms=[atom.copy() for atom in self.atoms], name=self.name)

    def residues(self) -> list[Residue]:
        order: list[tuple[str, int, str]] = []
        groups: dict[tuple[str, int, str], Residue] = {}
        for atom in self.atoms:
            key = (atom.chain, atom.resseq, atom.icode)
            residue = groups.get(key)
            if residue is None:
                residue = Residue(atom.resname, atom.chain, atom.resseq, atom.icode)
                groups[key] = residue
                order.append(key)
            residue.atoms[atom.name] = atom
        return [groups[key] for key in order]

    def ca_coordinates(self) -> np.ndarray:
        coords = []
        for residue in self.residues():
            atom = residue.atoms.get("CA")
            if atom is None:
                continue
            coords.append(atom.xyz)
        if not coords:
            return np.zeros((0, 3))
        return np.vstack(coords)

    def plddt_from_bfactors(self) -> np.ndarray | None:
        """AlphaFold writes the same pLDDT on every atom of a residue, in 0–100."""
        residues = self.residues()
        if not residues:
            return None
        values: list[float] = []
        constant = 0
        for residue in residues:
            factors = [atom.bfactor for atom in residue.atoms.values()]
            if not factors:
                return None
            if max(factors) - min(factors) <= 0.15:
                constant += 1
            ca = residue.atoms.get("CA")
            values.append(float(ca.bfactor if ca is not None else factors[0]))
        array = np.asarray(values, dtype=float)
        if constant < 0.9 * len(residues):
            return None
        if np.nanmin(array) < 0 or np.nanmax(array) > 100:
            return None
        return array


def parse_structure(path: str | Path) -> Structure:
    file_path = Path(path)
    text = file_path.read_text()
    if file_path.suffix.lower() in {".cif", ".mmcif"} or text.lstrip().startswith("data_"):
        atoms = _parse_cif(text)
    else:
        atoms = _parse_pdb(text)
    if not atoms:
        raise ValueError(f"No protein atoms found in {file_path}")
    return Structure(atoms=atoms, name=file_path.name)


def to_pdb(structure: Structure) -> str:
    lines = []
    for serial, atom in enumerate(structure.atoms, start=1):
        lines.append(_format_atom(serial, atom))
    lines.append("END")
    return "\n".join(lines) + "\n"


def _format_atom(serial: int, atom: Atom) -> str:
    name = _pdb_atom_name(atom.name)
    alt = (atom.altloc or " ")[:1]
    icode = (atom.icode or " ")[:1]
    element = (atom.element or element_from_name(atom.name))[:2]
    return (
        f"ATOM  {serial:5d} {name}{alt}{atom.resname:>3s} {atom.chain[:1] or 'A'}"
        f"{atom.resseq:4d}{icode}   "
        f"{atom.xyz[0]:8.3f}{atom.xyz[1]:8.3f}{atom.xyz[2]:8.3f}"
        f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}"
        f"          {element:>2s}"
    )


def _pdb_atom_name(name: str) -> str:
    name = name.strip().upper()
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3s}"


def _parse_pdb(text: str) -> list[Atom]:
    atoms: list[Atom] = []
    seen_model = False
    for raw in text.splitlines():
        record = raw[:6].strip().upper()
        if record == "MODEL":
            if seen_model and atoms:
                break
            seen_model = True
            continue
        if record == "ENDMDL" and atoms:
            break
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(raw) < 54:
            continue
        alt = raw[16:17].strip()
        if alt not in {"", "A"}:
            continue
        name = raw[12:16].strip().upper()
        resname = raw[17:20].strip().upper() or "UNK"
        if record == "HETATM" or resname in {"HOH", "WAT", "DOD"}:
            continue
        element = raw[76:78].strip().upper() if len(raw) >= 78 else ""
        if _is_hydrogen(name, element):
            continue
        try:
            xyz = np.array(
                [float(raw[30:38]), float(raw[38:46]), float(raw[46:54])],
                dtype=float,
            )
            occupancy = float(raw[54:60]) if raw[54:60].strip() else 1.0
            bfactor = float(raw[60:66]) if len(raw) >= 66 and raw[60:66].strip() else 0.0
            resseq = int(raw[22:26])
        except ValueError:
            continue
        atoms.append(
            Atom(
                name=name,
                element=element or element_from_name(name),
                resname=resname,
                chain=(raw[21:22].strip() or "A"),
                resseq=resseq,
                icode=raw[26:27].strip(),
                xyz=xyz,
                bfactor=bfactor,
                occupancy=occupancy,
                altloc=alt,
            )
        )
    return atoms


def _parse_cif(text: str) -> list[Atom]:
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            index += 1
            columns: list[str] = []
            while index < len(lines) and lines[index].strip().startswith("_"):
                columns.append(lines[index].strip())
                index += 1
            if not columns or not columns[0].startswith("_atom_site."):
                continue
            rows: list[list[str]] = []
            current: list[str] = []
            while index < len(lines):
                stripped = lines[index].strip()
                if not stripped or stripped.startswith("#"):
                    if current:
                        break
                    index += 1
                    continue
                if stripped.startswith(("loop_", "data_", "_")):
                    break
                if stripped.startswith(";"):
                    index += 1
                    blob = []
                    while index < len(lines) and not lines[index].startswith(";"):
                        blob.append(lines[index])
                        index += 1
                    current.append(" ".join(blob).strip())
                    index += 1
                else:
                    current.extend(_cif_tokens(stripped))
                    index += 1
                if len(current) >= len(columns):
                    rows.append(current[: len(columns)])
                    current = current[len(columns) :]
            return _atoms_from_cif_rows(columns, rows)
        index += 1
    return []


def _atoms_from_cif_rows(columns: list[str], rows: list[list[str]]) -> list[Atom]:
    def col(*names: str) -> int | None:
        for name in names:
            token = name if name.startswith("_") else f"_atom_site.{name}"
            if token in columns:
                return columns.index(token)
        return None

    group_i = col("group_PDB")
    atom_i = col("label_atom_id", "auth_atom_id")
    element_i = col("type_symbol")
    resname_i = col("label_comp_id", "auth_comp_id")
    chain_i = col("auth_asym_id", "label_asym_id")
    seq_i = col("auth_seq_id", "label_seq_id")
    icode_i = col("pdbx_PDB_ins_code")
    alt_i = col("label_alt_id")
    x_i = col("Cartn_x")
    y_i = col("Cartn_y")
    z_i = col("Cartn_z")
    occ_i = col("occupancy")
    b_i = col("B_iso_or_equiv")
    if None in (atom_i, x_i, y_i, z_i):
        return []

    atoms: list[Atom] = []
    for row in rows:
        if group_i is not None and row[group_i].upper() not in {"ATOM", "."}:
            continue
        name = row[atom_i].strip().upper()
        element = row[element_i].strip().upper() if element_i is not None else ""
        if element in {".", "?"}:
            element = ""
        if _is_hydrogen(name, element):
            continue
        resname = row[resname_i].strip().upper() if resname_i is not None else "UNK"
        if resname in {".", "?"}:
            resname = "UNK"
        if resname in {"HOH", "WAT", "DOD"}:
            continue
        alt = row[alt_i].strip() if alt_i is not None else ""
        if alt in {".", "?"}:
            alt = ""
        if alt not in {"", "A"}:
            continue
        chain = row[chain_i].strip() if chain_i is not None else "A"
        if chain in {"", ".", "?"}:
            chain = "A"
        icode = row[icode_i].strip() if icode_i is not None else ""
        if icode in {".", "?"}:
            icode = ""
        seq_token = row[seq_i].strip() if seq_i is not None else "1"
        if seq_token in {".", "?"}:
            continue
        try:
            xyz = np.array([float(row[x_i]), float(row[y_i]), float(row[z_i])], dtype=float)
            resseq = int(float(seq_token))
            occupancy = _cif_float(row[occ_i]) if occ_i is not None else 1.0
            bfactor = _cif_float(row[b_i]) if b_i is not None else 0.0
        except ValueError:
            continue
        atoms.append(
            Atom(
                name=name,
                element=element or element_from_name(name),
                resname=resname,
                chain=chain[:1] or "A",
                resseq=resseq,
                icode=icode,
                xyz=xyz,
                bfactor=0.0 if bfactor is None else bfactor,
                occupancy=1.0 if occupancy is None else occupancy,
                altloc=alt,
            )
        )
    return atoms


def _cif_float(token: str) -> float | None:
    if token in {".", "?"}:
        return None
    return float(token)


def _cif_tokens(line: str) -> list[str]:
    tokens: list[str] = []
    index = 0
    while index < len(line):
        if line[index].isspace():
            index += 1
            continue
        if line[index] in {"'", '"'}:
            quote = line[index]
            index += 1
            start = index
            while index < len(line) and line[index] != quote:
                index += 1
            tokens.append(line[start:index])
            index += 1
            continue
        start = index
        while index < len(line) and not line[index].isspace():
            index += 1
        tokens.append(line[start:index])
    return tokens


def _is_hydrogen(name: str, element: str) -> bool:
    if element in {"H", "D"}:
        return True
    return name.startswith("H")
