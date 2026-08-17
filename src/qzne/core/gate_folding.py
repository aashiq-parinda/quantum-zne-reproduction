"""Gate Folding for Noise Amplification — Core ZNE Technique.

The key physical insight of Temme et al. 2017 Section III.A:

Noise amplification via gate folding replaces each unitary U with an
equivalent operation that has the same ideal action but higher noise:

    U  →  U (U† U)^n    [amplifies noise by factor c = 2n+1]
    U  →  (U† U)^n U    [alternative folding order]

Since U†U = I exactly (unitaries are perfect), the ideal computation
is unchanged. But each additional gate application adds noise,
so the effective depolarizing strength becomes:

    λ_eff = c · λ    where  c = 2n+1 ∈ {1, 3, 5, 7, ...}

This gives us noise levels {λ, 3λ, 5λ, ...} by choosing n = 0, 1, 2, ...

For partial folding of a d-gate circuit:
    - Full folds: n complete U → UU†U replacements
    - Fractional: only first k gates folded

Reference: Temme, Bravyi, Gambetta (2017), PRL 119, 180509. Sec. III.A
           Giurgica-Tiron et al. (2020), arXiv:2005.10921 (gate folding formalization)
"""

import numpy as np
from typing import List, Tuple


def fold_gate(U: np.ndarray, n_folds: int = 1) -> np.ndarray:
    """Fold a unitary gate U → U (U† U)^n.

    Ideal action: Identical to U (since U†U = I).
    Noise: Amplified by factor (2n+1).

    Parameters
    ----------
    U : np.ndarray, shape (d, d) — unitary gate matrix
    n_folds : int
        Number of additional UU† pairs (n=0: identity fold, n=1: U→UU†U)

    Returns
    -------
    np.ndarray — effective unitary (mathematically U, physically c·λ noisy)
    """
    if n_folds < 0:
        raise ValueError(f"n_folds={n_folds} must be ≥ 0")
    result = U.copy()
    for _ in range(n_folds):
        result = result @ U.conj().T @ U
    return result


def noise_scale_factor(n_folds: int) -> int:
    """Return noise scale factor c = 2n+1 for n gate fold repetitions."""
    return 2 * n_folds + 1


def build_noise_levels(base_lambda: float, n_folds_list: List[int]) -> List[float]:
    """Build list of amplified noise levels from base λ and fold counts.

    λ_i = c_i · base_lambda    where c_i = 2n_i + 1

    Parameters
    ----------
    base_lambda : float — physical device noise rate
    n_folds_list : list of int — fold counts e.g. [0, 1, 2] → scales [1, 3, 5]

    Returns
    -------
    list of effective noise rates [λ_1, λ_2, ..., λ_k]
    """
    return [noise_scale_factor(n) * base_lambda for n in n_folds_list]


def global_fold_circuit(
    gates: List[np.ndarray],
    n_folds: int,
) -> List[np.ndarray]:
    """Apply global folding to an entire circuit: G → G (G† G)^n.

    For full circuit G = U_d ... U_1:
    G → G · (G† G)^n = G · (U_1† ... U_d† · U_d ... U_1)^n

    Practically implemented as reversing + repeating the gate list.

    Parameters
    ----------
    gates : list of gate matrices [U_1, U_2, ..., U_d] (ordered)
    n_folds : int

    Returns
    -------
    list of gates with length d * (2n+1)
    """
    folded = list(gates)
    for _ in range(n_folds):
        folded += [g.conj().T for g in reversed(gates)]
        folded += list(gates)
    return folded


def local_fold_circuit(
    gates: List[np.ndarray],
    n_folds: int,
    k: int | None = None,
) -> List[np.ndarray]:
    """Apply local gate-level folding, optionally folding only first k gates.

    Each gate U_i → U_i (U_i† U_i)^n.

    Parameters
    ----------
    gates : list of gate matrices
    n_folds : int — fold repetitions per gate
    k : int or None — if set, only fold first k gates (partial fold)

    Returns
    -------
    list of folded gates
    """
    k = k if k is not None else len(gates)
    result = []
    for i, g in enumerate(gates):
        if i < k:
            result.append(fold_gate(g, n_folds))
        else:
            result.append(g)
    return result


def fractional_scale_circuit(
    gates: List[np.ndarray],
    target_scale: float,
) -> List[np.ndarray]:
    """Produce a circuit with noise scaled by an arbitrary factor s ≥ 1.

    For non-integer scale s between consecutive integer scales 2n+1 and 2n+3:
      - Apply n full circuit folds globally
      - Apply 1 additional fold to the first k gates, where k = round((s - (2n+1)) * d / 2)

    Parameters
    ----------
    gates : list of gate matrices
    target_scale : float — desired noise scale factor ≥ 1.0

    Returns
    -------
    list of gates achieving approximately target_scale × base noise
    """
    if target_scale < 1.0:
        raise ValueError(f"target_scale={target_scale} must be ≥ 1.0")
    d = len(gates)
    n_full = int((target_scale - 1) // 2)
    remainder = target_scale - (2 * n_full + 1)
    k_partial = int(round(remainder * d / 2))

    folded = global_fold_circuit(gates, n_full)
    if k_partial > 0:
        folded = local_fold_circuit(folded, 1, k=k_partial)
    return folded
