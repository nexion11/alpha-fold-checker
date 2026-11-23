"""Read AlphaFold pLDDT and PAE JSON."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def load_confidence_file(path: str | Path) -> dict:
    data = json.loads(Path(path).read_text())
    plddt = _find_plddt(data)
    pae = _find_pae(data)
    ptm = _find_ptm(data)
    if plddt is None and pae is None and ptm is None:
        raise ValueError(f"No pLDDT, PAE, or pTM value found in {path}")
    return {"plddt": plddt, "pae": pae, "ptm": ptm}


def _find_plddt(data) -> np.ndarray | None:
    if isinstance(data, dict):
        for key in ("plddt", "pLDDT"):
            if key in data:
                return _as_vector(data[key])
        for value in data.values():
            found = _find_plddt(value)
            if found is not None:
                return found
    if isinstance(data, list) and data and isinstance(data[0], dict):
        return _find_plddt(data[0])
    return None


def _find_ptm(data) -> float | None:
    if isinstance(data, dict):
        for key in ("ptm", "pTM", "iptm"):
            value = data.get(key)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                return float(value)
    if isinstance(data, list) and data and isinstance(data[0], dict):
        return _find_ptm(data[0])
    return None


def _find_pae(data) -> np.ndarray | None:
    if isinstance(data, dict):
        for key in ("predicted_aligned_error", "pae", "PAE"):
            if key in data:
                matrix = _as_matrix(data[key])
                if matrix is not None:
                    return matrix
    if isinstance(data, list) and data and isinstance(data[0], dict):
        record = data[0]
        if "distance" in record and "residue1" in record and "residue2" in record:
            residue1 = np.asarray(record["residue1"], dtype=int)
            residue2 = np.asarray(record["residue2"], dtype=int)
            distance = np.asarray(record["distance"], dtype=float)
            size = int(max(residue1.max(), residue2.max()))
            matrix = np.full((size, size), np.nan)
            matrix[residue1 - 1, residue2 - 1] = distance
            return matrix
        return _find_pae(data[0])
    return _as_matrix(data)


def _as_vector(value) -> np.ndarray | None:
    if isinstance(value, list) and value and isinstance(value[0], list):
        value = value[0]
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or array.size < 2:
        return None
    return array


def _as_matrix(value) -> np.ndarray | None:
    if isinstance(value, list) and value and isinstance(value[0], dict):
        return None
    array = np.asarray(value, dtype=float)
    if array.ndim == 3:
        array = array[0]
    if array.ndim != 2 or array.shape[0] != array.shape[1] or array.shape[0] < 2:
        return None
    return array
