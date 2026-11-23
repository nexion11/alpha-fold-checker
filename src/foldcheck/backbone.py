"""Build ideal peptide backbones from phi/psi angles."""

from __future__ import annotations

import numpy as np

from foldcheck.geometry import (
    ANGLE_IDEAL,
    BOND_IDEAL,
    PEPTIDE_ANGLE_C_N_CA,
    PEPTIDE_ANGLE_CA_C_N,
    PEPTIDE_BOND,
)
from foldcheck.structure import Atom, Structure

HELIX = (-57.0, -47.0)
STRAND = (-119.0, 113.0)
CB_DIHEDRAL = 122.6
CB_ANGLE = 110.4


def build_peptide(
    length: int,
    phi: float = HELIX[0],
    psi: float = HELIX[1],
    chain: str = "A",
    resname: str = "ALA",
) -> Structure:
    if length < 2:
        raise ValueError("A peptide needs at least 2 residues")
    prev_n = np.array([0.0, 0.0, 0.0])
    prev_ca = np.array([BOND_IDEAL[("N", "CA")], 0.0, 0.0])
    dummy = np.array([0.0, 1.0, 0.0])
    prev_c = place_atom(
        dummy,
        prev_n,
        prev_ca,
        BOND_IDEAL[("CA", "C")],
        ANGLE_IDEAL[("N", "CA", "C")],
        phi,
    )
    atoms: list[Atom] = []
    _append_residue(atoms, 1, chain, resname, prev_n, prev_ca, prev_c)
    for resseq in range(2, length + 1):
        prev_n = place_atom(prev_n, prev_ca, prev_c, PEPTIDE_BOND, PEPTIDE_ANGLE_CA_C_N, psi)
        prev_ca = place_atom(
            prev_ca,
            prev_c,
            prev_n,
            BOND_IDEAL[("N", "CA")],
            PEPTIDE_ANGLE_C_N_CA,
            180.0,
        )
        prev_c = place_atom(
            prev_c,
            prev_n,
            prev_ca,
            BOND_IDEAL[("CA", "C")],
            ANGLE_IDEAL[("N", "CA", "C")],
            phi,
        )
        _append_residue(atoms, resseq, chain, resname, prev_n, prev_ca, prev_c)
    return Structure(atoms=atoms, name=f"{resname}{length}")


def place_atom(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    length: float,
    angle_deg: float,
    dihedral_deg: float,
) -> np.ndarray:
    """Place atom D given A, B, C, the C–D length, the B–C–D angle, and dihedral A–B–C–D."""
    angle = np.deg2rad(angle_deg)
    dihedral = np.deg2rad(dihedral_deg)
    bc = c - b
    bc_hat = bc / np.linalg.norm(bc)
    normal = np.cross(b - a, bc)
    normal_length = np.linalg.norm(normal)
    if normal_length < 1e-8:
        raise ValueError("Cannot place an atom from three colinear points")
    normal_hat = normal / normal_length
    in_plane = np.cross(normal_hat, bc_hat)
    return (
        c
        - bc_hat * length * np.cos(angle)
        + in_plane * length * np.sin(angle) * np.cos(dihedral)
        + normal_hat * length * np.sin(angle) * np.sin(dihedral)
    )


def _append_residue(
    atoms: list[Atom],
    resseq: int,
    chain: str,
    resname: str,
    nitrogen: np.ndarray,
    ca: np.ndarray,
    carbon: np.ndarray,
) -> None:
    oxygen = place_atom(nitrogen, ca, carbon, BOND_IDEAL[("C", "O")], ANGLE_IDEAL[("CA", "C", "O")], 180.0)
    cb = place_atom(carbon, nitrogen, ca, BOND_IDEAL[("CA", "CB")], CB_ANGLE, CB_DIHEDRAL)
    records = (
        ("N", "N", nitrogen),
        ("CA", "C", ca),
        ("C", "C", carbon),
        ("O", "O", oxygen),
        ("CB", "C", cb),
    )
    for name, element, xyz in records:
        atoms.append(
            Atom(
                name=name,
                element=element,
                resname=resname,
                chain=chain,
                resseq=resseq,
                icode="",
                xyz=np.array(xyz, dtype=float, copy=True),
                bfactor=0.0,
                occupancy=1.0,
            )
        )
