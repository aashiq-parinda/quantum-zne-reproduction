"""Noise channel models for ZNE reproduction.

Implements Kraus operator representations for:
  - Depolarizing channel (used in Temme et al. 2017)
  - Amplitude damping / T₁ relaxation
  - Dephasing / T₂ channel

Depolarizing channel at strength λ:
    ε_λ(ρ) = (1 - λ)ρ + (λ/3)(XρX + YρY + ZρZ)
           = (1 - 4λ/3)ρ + (4λ/3) I/2

Kraus operators:
    K₀ = √(1 - λ) I
    K₁ = √(λ/3) X
    K₂ = √(λ/3) Y
    K₃ = √(λ/3) Z
"""

import numpy as np


# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def depolarizing_kraus(lam: float) -> list[np.ndarray]:
    """Return single-qubit depolarizing Kraus operators at noise rate λ.

    ε_λ(ρ) = (1-λ)ρ + (λ/3)(XρX + YρY + ZρZ)

    Parameterization matches Temme et al. 2017 convention.

    Parameters
    ----------
    lam : float
        Depolarizing strength λ ∈ [0, 3/4].

    Returns
    -------
    list of 4 Kraus operators K₀..K₃
    """
    if lam < 0 or lam > 0.75:
        raise ValueError(f"λ={lam} must be in [0, 0.75] for valid depolarizing channel")
    K0 = np.sqrt(1.0 - lam) * I2
    K1 = np.sqrt(lam / 3.0) * X
    K2 = np.sqrt(lam / 3.0) * Y
    K3 = np.sqrt(lam / 3.0) * Z
    return [K0, K1, K2, K3]


def dephasing_kraus(lam: float) -> list[np.ndarray]:
    """Return single-qubit dephasing (phase-flip) Kraus operators.

    ε_λ(ρ) = (1-λ)ρ + λ ZρZ

    Parameters
    ----------
    lam : float
        Dephasing strength λ ∈ [0, 0.5].
    """
    if lam < 0 or lam > 0.5:
        raise ValueError(f"λ={lam} must be in [0, 0.5] for dephasing channel")
    K0 = np.sqrt(1.0 - lam) * I2
    K1 = np.sqrt(lam) * Z
    return [K0, K1]


def amplitude_damping_kraus(gamma: float) -> list[np.ndarray]:
    """Return single-qubit amplitude damping Kraus operators (T₁ relaxation).

    K₀ = [[1, 0], [0, √(1-γ)]]
    K₁ = [[0, √γ], [0, 0]]

    Parameters
    ----------
    gamma : float
        Damping parameter γ ∈ [0, 1] related to T₁ via γ = 1 - exp(-t/T₁).
    """
    if gamma < 0 or gamma > 1:
        raise ValueError(f"γ={gamma} must be in [0, 1]")
    K0 = np.array([[1, 0], [0, np.sqrt(1.0 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]


def apply_kraus_channel(
    rho: np.ndarray,
    kraus_ops: list[np.ndarray],
) -> np.ndarray:
    """Apply a Kraus channel ε to density matrix ρ.

    ε(ρ) = Σ_k K_k ρ K_k†

    Parameters
    ----------
    rho : np.ndarray, shape (2, 2) — single-qubit density matrix
    kraus_ops : list of Kraus operators

    Returns
    -------
    np.ndarray — noisy density matrix ε(ρ)
    """
    result = np.zeros_like(rho)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return result


def verify_completeness(kraus_ops: list[np.ndarray], tol: float = 1e-10) -> bool:
    """Verify Σ_k K_k† K_k = I (completeness relation for trace-preserving map)."""
    n = kraus_ops[0].shape[0]
    check = sum(K.conj().T @ K for K in kraus_ops)
    return np.allclose(check, np.eye(n), atol=tol)
