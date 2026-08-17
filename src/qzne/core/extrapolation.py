"""Richardson Extrapolation & Polynomial Fitting for ZNE.

The mathematical engine of Zero Noise Extrapolation:

Given measurements ⟨O⟩_{λ_i} at noise levels λ_1 < λ_2 < ... < λ_k,
extrapolate to the zero-noise limit ⟨O⟩_0.

METHOD 1 — Linear Extrapolation (2-point Richardson):
    Assumes ⟨O⟩_λ ≈ ⟨O⟩_0 + a₁λ (linear noise model)
    From two points (λ₁, E₁) and (λ₂, E₂):
        ⟨O⟩_0 = (λ₂·E₁ - λ₁·E₂) / (λ₂ - λ₁)

METHOD 2 — Richardson Extrapolation (n-point, order n-1 error cancellation):
    Barycentric form of Richardson's deferred approach:
        ⟨O⟩_0 = Σ_k c_k · ⟨O⟩_{λ_k}
    where   c_k = Π_{j≠k} λ_j / (λ_j - λ_k)     (Richardson coefficients)

    For equally-spaced scales c_i = (2i-1) with i=1..n (standard folding):
        λ_i = (2i-1)λ,   so c_k = Π_{j≠k} (2j-1) / ((2j-1)-(2k-1))

METHOD 3 — Polynomial Fitting (degree m):
    Fit ⟨O⟩_λ = a₀ + a₁λ + ... + a_mλ^m to all data points via least squares.
    Extrapolated value = a₀.

Reference: Temme, Bravyi, Gambetta (2017), PRL 119, 180509. Eq. (8)-(10)
"""

import numpy as np
from typing import Tuple


def richardson_coefficients(lambdas: np.ndarray) -> np.ndarray:
    """Compute Richardson extrapolation coefficients c_k for given noise levels.

    c_k = Π_{j≠k} λ_j / (λ_j - λ_k)

    These coefficients are such that Σ_k c_k · ⟨O⟩_{λ_k} = ⟨O⟩_0
    to order n-1 in λ, where n = len(lambdas).

    Parameters
    ----------
    lambdas : np.ndarray — noise levels [λ₁, λ₂, ..., λ_n]

    Returns
    -------
    np.ndarray — Richardson coefficients [c₁, ..., c_n]
    """
    n = len(lambdas)
    coeffs = np.ones(n)
    for k in range(n):
        for j in range(n):
            if j != k:
                denom = lambdas[j] - lambdas[k]
                if abs(denom) < 1e-14:
                    raise ValueError(f"Duplicate noise levels at i={j},{k}: λ={lambdas[j]}")
                coeffs[k] *= lambdas[j] / denom
    return coeffs


def richardson_extrapolate(
    lambdas: np.ndarray,
    expectations: np.ndarray,
) -> float:
    """Zero-noise extrapolation via Richardson's deferred approach.

    ⟨O⟩_0 ≈ Σ_k c_k · ⟨O⟩_{λ_k}

    Parameters
    ----------
    lambdas : np.ndarray — noise levels
    expectations : np.ndarray — measured expectation values at each λ

    Returns
    -------
    float — Richardson-extrapolated zero-noise expectation value
    """
    coeffs = richardson_coefficients(lambdas)
    return float(np.dot(coeffs, expectations))


def linear_extrapolate(
    lam1: float, exp1: float,
    lam2: float, exp2: float,
) -> float:
    """Linear 2-point zero-noise extrapolation.

    ⟨O⟩_0 = (λ₂·E₁ - λ₁·E₂) / (λ₂ - λ₁)

    Parameters
    ----------
    lam1, exp1 : first noise level and expectation
    lam2, exp2 : second noise level and expectation
    """
    return (lam2 * exp1 - lam1 * exp2) / (lam2 - lam1)


def polynomial_extrapolate(
    lambdas: np.ndarray,
    expectations: np.ndarray,
    degree: int | None = None,
) -> Tuple[float, np.ndarray]:
    """Polynomial least-squares extrapolation to λ=0.

    Fits ⟨O⟩_λ = a₀ + a₁λ + ... + a_m λ^m via numpy polyfit.

    Parameters
    ----------
    lambdas : np.ndarray — noise levels
    expectations : np.ndarray — expectation values at each level
    degree : int or None
        Polynomial degree. If None, uses len(lambdas)-1 (interpolation).

    Returns
    -------
    (extrapolated_value, coefficients) — a₀ = extrapolated zero-noise value
    """
    if degree is None:
        degree = len(lambdas) - 1
    degree = min(degree, len(lambdas) - 1)

    coeffs = np.polyfit(lambdas, expectations, degree)
    poly = np.poly1d(coeffs)
    return float(poly(0.0)), coeffs


def exponential_extrapolate(
    lambdas: np.ndarray,
    expectations: np.ndarray,
) -> float:
    """Exponential model extrapolation: ⟨O⟩_λ = A·exp(-B·λ) + C.

    Fits the log-transformed data to extract A, B, C and extrapolate to λ=0.
    More physically motivated for amplitude-damping noise channels.

    Parameters
    ----------
    lambdas : np.ndarray — noise levels
    expectations : np.ndarray — expectation values at each level

    Returns
    -------
    float — A + C (extrapolated zero-noise expectation)
    """
    # Shift to make all values positive for log transform
    min_val = expectations.min()
    offset = abs(min_val) + 0.01 if min_val <= 0 else 0.0
    shifted = expectations + offset

    # Fit log(E + offset) = log(A) - B*λ
    log_shifted = np.log(shifted)
    coeffs = np.polyfit(lambdas, log_shifted, 1)
    B = -coeffs[0]
    log_A = coeffs[1]
    A = np.exp(log_A)

    return float(A - offset)
