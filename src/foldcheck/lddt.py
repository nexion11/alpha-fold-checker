"""C-alpha lDDT and RMSD.

lDDT compares pairwise distances, so it does not depend on superposition.
It is the same family of score that AlphaFold's pLDDT is trained to predict.
"""

from __future__ import annotations

import numpy as np

LDDT_CUTOFF = 15.0
LDDT_THRESHOLDS = (0.5, 1.0, 2.0, 4.0)


def lddt_ca(model: np.ndarray, reference: np.ndarray, cutoff: float = LDDT_CUTOFF) -> float:
    model = np.asarray(model, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if model.shape != reference.shape or model.ndim != 2 or model.shape[1] != 3:
        raise ValueError("model and reference CA coordinates must share shape (n, 3)")
    count = model.shape[0]
    if count < 2:
        return float("nan")
    scores: list[float] = []
    for start in range(0, count, 256):
        stop = min(count, start + 256)
        ref_dist = np.linalg.norm(reference[start:stop, None, :] - reference[None, :, :], axis=-1)
        model_dist = np.linalg.norm(model[start:stop, None, :] - model[None, :, :], axis=-1)
        neighbors = ref_dist < cutoff
        for row in range(stop - start):
            neighbors[row, start + row] = False
        included = neighbors.sum(axis=1)
        delta = np.abs(model_dist - ref_dist)
        for row in range(stop - start):
            if included[row] == 0:
                continue
            row_delta = delta[row, neighbors[row]]
            preserved = [(row_delta < threshold).mean() for threshold in LDDT_THRESHOLDS]
            scores.append(float(np.mean(preserved)))
    if not scores:
        return float("nan")
    return float(np.mean(scores))


def kabsch_rmsd(model: np.ndarray, reference: np.ndarray) -> float:
    model = np.asarray(model, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if model.shape != reference.shape or model.shape[0] < 3:
        raise ValueError("RMSD needs matching coordinates and at least 3 atoms")
    centered_model = model - model.mean(axis=0)
    centered_reference = reference - reference.mean(axis=0)
    covariance = centered_model.T @ centered_reference
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    v_matrix = vt_matrix.T
    sign = np.linalg.det(v_matrix @ u_matrix.T)
    correction = np.diag([1.0, 1.0, float(sign)])
    rotation = v_matrix @ correction @ u_matrix.T
    aligned = centered_model @ rotation
    delta = aligned - centered_reference
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
