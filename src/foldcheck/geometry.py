"""Backbone geometry checks for a predicted structure.

Thresholds are simpler than MolProbity. They catch broken covalent geometry,
hard steric clashes, chain breaks, and clearly disallowed phi/psi angles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from foldcheck.structure import Atom, Residue, Structure

BOND_IDEAL = {
    ("N", "CA"): 1.458,
    ("CA", "C"): 1.525,
    ("C", "O"): 1.231,
    ("CA", "CB"): 1.530,
}
PEPTIDE_BOND = 1.329
BOND_TOLERANCE = 0.12

ANGLE_IDEAL = {
    ("N", "CA", "C"): 111.2,
    ("CA", "C", "O"): 120.8,
}
PEPTIDE_ANGLE_CA_C_N = 116.2
PEPTIDE_ANGLE_C_N_CA = 121.7
ANGLE_TOLERANCE = 10.0

CA_DISTANCE = 3.80
HEAVY_CLASH = 2.20
CA_ONLY_CLASH = 3.00
CHAIN_BREAK_MIN = 2.90
CHAIN_BREAK_MAX = 4.50


@dataclass
class GeometryReport:
    clashes: int
    bond_outlier_fraction: float | None
    angle_outlier_fraction: float | None
    rama_outlier_fraction: float | None
    chain_breaks: int
    ca_dist_deviation: float | None
    radius_of_gyration: float | None
    residue_count: int
    chain_count: int

    def score(self) -> float:
        residue_count = max(self.residue_count, 1)
        penalty = min((self.clashes / residue_count) / 0.15, 1.0) * 0.50
        if self.bond_outlier_fraction is not None:
            penalty += min(self.bond_outlier_fraction / 0.10, 1.0) * 0.20
        if self.angle_outlier_fraction is not None:
            penalty += min(self.angle_outlier_fraction / 0.15, 1.0) * 0.10
        if self.rama_outlier_fraction is not None:
            penalty += min(self.rama_outlier_fraction / 0.30, 1.0) * 0.15
        penalty += min(self.chain_breaks * 0.08, 0.25)
        return float(np.clip(1.0 - penalty, 0.0, 1.0))


def element_from_name(name: str) -> str:
    letters = "".join(character for character in name.strip().upper() if character.isalpha())
    if not letters:
        return "C"
    if letters.startswith(("CA", "CB", "CG", "CD", "CE", "CZ", "CH")):
        return "C"
    return letters[0]


def analyze_geometry(structure: Structure) -> GeometryReport:
    residues = structure.residues()
    chains = {residue.chain for residue in residues}
    clashes = _count_clashes(residues)
    bond_fraction = _fraction(_bond_checks(residues))
    angle_fraction = _fraction(_angle_checks(residues))
    rama_fraction = _fraction(_rama_checks(residues))
    breaks, ca_deviation = _ca_trace(residues)
    return GeometryReport(
        clashes=clashes,
        bond_outlier_fraction=bond_fraction,
        angle_outlier_fraction=angle_fraction,
        rama_outlier_fraction=rama_fraction,
        chain_breaks=breaks,
        ca_dist_deviation=ca_deviation,
        radius_of_gyration=_radius_of_gyration(residues),
        residue_count=len(residues),
        chain_count=len(chains),
    )


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def rama_favored(phi: float, psi: float, resname: str) -> bool:
    if resname.upper() == "GLY":
        return True
    if -180 <= phi <= -40 and (80 <= psi <= 180 or -180 <= psi <= -160):
        return True
    if -160 <= phi <= -40 and -80 <= psi <= 50:
        return True
    if 40 <= phi <= 90 and -20 <= psi <= 80:
        return True
    return False


def _fraction(checks: list[bool]) -> float | None:
    if not checks:
        return None
    return float(sum(1 for flag in checks if flag) / len(checks))


def _bond_checks(residues: list[Residue]) -> list[bool]:
    outliers: list[bool] = []
    for residue in residues:
        for (left, right), ideal in BOND_IDEAL.items():
            pair = _pair(residue, residue, left, right)
            if pair is None:
                continue
            outliers.append(abs(pair - ideal) > BOND_TOLERANCE)
    for first, second in _consecutive(residues):
        distance = _pair(first, second, "C", "N")
        if distance is None:
            continue
        outliers.append(abs(distance - PEPTIDE_BOND) > BOND_TOLERANCE)
    return outliers


def _angle_checks(residues: list[Residue]) -> list[bool]:
    outliers: list[bool] = []
    for residue in residues:
        for (a_name, b_name, c_name), ideal in ANGLE_IDEAL.items():
            angle = _angle(residue, residue, residue, a_name, b_name, c_name)
            if angle is None:
                continue
            outliers.append(abs(angle - ideal) > ANGLE_TOLERANCE)
    for first, second in _consecutive(residues):
        peptide = _angle(first, first, second, "CA", "C", "N")
        if peptide is not None:
            outliers.append(abs(peptide - PEPTIDE_ANGLE_CA_C_N) > ANGLE_TOLERANCE)
        incoming = _angle(first, second, second, "C", "N", "CA")
        if incoming is not None:
            outliers.append(abs(incoming - PEPTIDE_ANGLE_C_N_CA) > ANGLE_TOLERANCE)
    return outliers


def _rama_checks(residues: list[Residue]) -> list[bool]:
    outliers: list[bool] = []
    for index, residue in enumerate(residues):
        if index == 0 or index + 1 >= len(residues):
            continue
        previous, following = residues[index - 1], residues[index + 1]
        if not _bonded(previous, residue) or not _bonded(residue, following):
            continue
        needed = ("C", "N", "CA", "C")
        atoms = (
            previous.atoms.get("C"),
            residue.atoms.get("N"),
            residue.atoms.get("CA"),
            residue.atoms.get("C"),
            following.atoms.get("N"),
        )
        if any(atom is None for atom in atoms):
            continue
        prev_c, nitrogen, ca, carbon, next_n = atoms
        phi = dihedral(prev_c.xyz, nitrogen.xyz, ca.xyz, carbon.xyz)
        psi = dihedral(nitrogen.xyz, ca.xyz, carbon.xyz, next_n.xyz)
        outliers.append(not rama_favored(phi, psi, residue.resname))
    return outliers


def _ca_trace(residues: list[Residue]) -> tuple[int, float | None]:
    breaks = 0
    deviations: list[float] = []
    for first, second in _consecutive(residues):
        left = first.atoms.get("CA")
        right = second.atoms.get("CA")
        if left is None or right is None:
            continue
        distance = float(np.linalg.norm(left.xyz - right.xyz))
        if distance < CHAIN_BREAK_MIN or distance > CHAIN_BREAK_MAX:
            breaks += 1
            continue
        deviations.append(abs(distance - CA_DISTANCE))
    if not deviations:
        return breaks, None
    return breaks, float(np.mean(deviations))


def _radius_of_gyration(residues: list[Residue]) -> float | None:
    coords = [residue.atoms["CA"].xyz for residue in residues if "CA" in residue.atoms]
    if len(coords) < 2:
        return None
    points = np.vstack(coords)
    centered = points - points.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))


def _count_clashes(residues: list[Residue]) -> int:
    atoms: list[Atom] = []
    residue_index: list[int] = []
    heavy_per_residue = []
    for index, residue in enumerate(residues):
        heavy = [
            atom
            for atom in residue.atoms.values()
            if atom.element not in {"H", "D"} and not atom.name.startswith("H")
        ]
        heavy_per_residue.append(len(heavy))
        for atom in heavy:
            atoms.append(atom)
            residue_index.append(index)
    if len(atoms) < 2:
        return 0
    ca_only = sum(1 for count in heavy_per_residue if count <= 1) >= 0.8 * len(residues)
    threshold = CA_ONLY_CLASH if ca_only else HEAVY_CLASH
    coords = np.vstack([atom.xyz for atom in atoms])
    owners = np.asarray(residue_index, dtype=int)
    return _grid_clashes(coords, owners, threshold)


def _grid_clashes(coords: np.ndarray, owners: np.ndarray, threshold: float) -> int:
    cells = np.floor(coords / threshold).astype(np.int64)
    buckets: dict[tuple[int, int, int], list[int]] = {}
    for index, cell in enumerate(cells):
        buckets.setdefault((int(cell[0]), int(cell[1]), int(cell[2])), []).append(index)
    clashes = 0
    limit = threshold * threshold
    for (x_cell, y_cell, z_cell), members in buckets.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbors = buckets.get((x_cell + dx, y_cell + dy, z_cell + dz))
                    if not neighbors:
                        continue
                    for left in members:
                        for right in neighbors:
                            if right <= left:
                                continue
                            if abs(int(owners[left]) - int(owners[right])) < 2:
                                continue
                            delta = coords[left] - coords[right]
                            if float(delta @ delta) < limit:
                                clashes += 1
    return clashes


def _consecutive(residues: list[Residue]) -> list[tuple[Residue, Residue]]:
    pairs = []
    for first, second in zip(residues, residues[1:]):
        if _bonded(first, second):
            pairs.append((first, second))
    return pairs


def _bonded(first: Residue, second: Residue) -> bool:
    return (
        first.chain == second.chain
        and first.icode == second.icode
        and second.resseq - first.resseq == 1
    )


def _pair(left: Residue, right: Residue, left_name: str, right_name: str) -> float | None:
    a = left.atoms.get(left_name)
    b = right.atoms.get(right_name)
    if a is None or b is None:
        return None
    return float(np.linalg.norm(a.xyz - b.xyz))


def _angle(
    a_res: Residue,
    b_res: Residue,
    c_res: Residue,
    a_name: str,
    b_name: str,
    c_name: str,
) -> float | None:
    a = a_res.atoms.get(a_name)
    b = b_res.atoms.get(b_name)
    c = c_res.atoms.get(c_name)
    if a is None or b is None or c is None:
        return None
    u = a.xyz - b.xyz
    v = c.xyz - b.xyz
    denom = np.linalg.norm(u) * np.linalg.norm(v)
    if denom == 0:
        return None
    cosine = float(np.clip(np.dot(u, v) / denom, -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))
