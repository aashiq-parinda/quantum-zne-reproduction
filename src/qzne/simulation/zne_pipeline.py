"""Full Zero Noise Extrapolation Pipeline.

Orchestrates the complete ZNE workflow from Temme et al. 2017:

1. Execute circuit at noise level λ      → ⟨O⟩_λ
2. Fold circuit to amplify noise: 3λ, 5λ → ⟨O⟩_{3λ}, ⟨O⟩_{5λ}
3. Apply Richardson / polynomial extrapolation → ⟨O⟩_0 (estimated)

Compare with ideal (λ=0) simulation to compute mitigation factor.
"""

import numpy as np
from typing import Callable, Dict, Any, List, Tuple
from qzne.core.gate_folding import noise_scale_factor
from qzne.core.extrapolation import (
    richardson_extrapolate,
    polynomial_extrapolate,
    linear_extrapolate,
    exponential_extrapolate,
)
from qzne.simulation.noisy_simulator import DensityMatrixSimulator, build_zz_observable


def run_circuit_at_noise(
    circuit_fn: Callable[[DensityMatrixSimulator, int], None],
    n_qubits: int,
    noise_lambda: float,
    noise_scale: int,
    observable: np.ndarray,
) -> float:
    """Run a circuit function at amplified noise level.

    Gate folding amplifies physical noise λ by factor c = 2n+1:
    Effective noise: λ_eff = c × λ

    Parameters
    ----------
    circuit_fn : callable(sim, n_folds) → None
        Function that applies the circuit to a DensityMatrixSimulator.
        Should apply each gate n_fold times extra.
    n_qubits : int
    noise_lambda : float — physical device noise rate
    noise_scale : int — noise amplification factor c (e.g. 1, 3, 5)
    observable : np.ndarray — Hermitian operator to measure

    Returns
    -------
    float — ⟨O⟩_{c·λ}
    """
    n_folds = (noise_scale - 1) // 2  # c = 2n+1 → n = (c-1)/2
    sim = DensityMatrixSimulator(
        n_qubits=n_qubits,
        noise_lambda=noise_lambda * noise_scale,
    )
    circuit_fn(sim, n_folds)
    return sim.expectation_value(observable)


def zne_workflow(
    circuit_fn: Callable[[DensityMatrixSimulator, int], None],
    n_qubits: int,
    base_lambda: float,
    noise_scales: List[int],
    observable: np.ndarray,
    method: str = "richardson",
) -> Dict[str, Any]:
    """Execute full ZNE workflow.

    Parameters
    ----------
    circuit_fn : callable — circuit builder
    n_qubits : int
    base_lambda : float — physical device noise rate
    noise_scales : list of int — noise amplification factors c [1, 3, 5, ...]
    observable : np.ndarray — Hermitian operator
    method : str — 'richardson', 'linear', 'polynomial', 'exponential'

    Returns
    -------
    dict with:
        - lambdas: noise levels used
        - expectations: ⟨O⟩ at each noise level
        - zne_estimate: extrapolated ⟨O⟩_0
        - ideal_value: noise-free reference
        - raw_error: |⟨O⟩_{λ} - ⟨O⟩_0|
        - zne_error: |⟨O⟩_ZNE - ⟨O⟩_0|
        - mitigation_factor: raw_error / zne_error
    """
    # Ideal (λ=0) reference
    sim_ideal = DensityMatrixSimulator(n_qubits=n_qubits, noise_lambda=0.0)
    circuit_fn(sim_ideal, 0)
    ideal_value = sim_ideal.expectation_value(observable)

    # Noisy executions at each amplified level
    lambdas = np.array([base_lambda * c for c in noise_scales])
    expectations = np.array([
        run_circuit_at_noise(circuit_fn, n_qubits, base_lambda, c, observable)
        for c in noise_scales
    ])

    # Extrapolation
    if method == "richardson":
        zne_estimate = richardson_extrapolate(lambdas, expectations)
    elif method == "linear":
        zne_estimate = linear_extrapolate(lambdas[0], expectations[0], lambdas[1], expectations[1])
    elif method == "polynomial":
        zne_estimate, _ = polynomial_extrapolate(lambdas, expectations)
    elif method == "exponential":
        zne_estimate = exponential_extrapolate(lambdas, expectations)
    else:
        raise ValueError(f"Unknown method '{method}'. Use: 'richardson', 'linear', 'polynomial', 'exponential'")

    raw_error = abs(expectations[0] - ideal_value)
    zne_error = abs(zne_estimate - ideal_value)
    mitigation_factor = raw_error / max(zne_error, 1e-14)

    return {
        "lambdas": lambdas,
        "expectations": expectations,
        "ideal_value": ideal_value,
        "zne_estimate": zne_estimate,
        "raw_error": raw_error,
        "zne_error": zne_error,
        "mitigation_factor": mitigation_factor,
        "method": method,
    }


def scan_noise_levels(
    circuit_fn: Callable,
    n_qubits: int,
    observable: np.ndarray,
    lambdas: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """Scan ⟨O⟩ vs λ for figure reproduction (replicates Fig 1 of Temme 2017).

    Parameters
    ----------
    circuit_fn : callable
    n_qubits : int
    observable : np.ndarray
    lambdas : np.ndarray — array of noise levels

    Returns
    -------
    (expectations, ideal_value) — arrays ready for plotting
    """
    expectations = []
    for lam in lambdas:
        sim = DensityMatrixSimulator(n_qubits=n_qubits, noise_lambda=lam)
        circuit_fn(sim, 0)
        expectations.append(sim.expectation_value(observable))

    sim_ideal = DensityMatrixSimulator(n_qubits=n_qubits, noise_lambda=0.0)
    circuit_fn(sim_ideal, 0)
    ideal = sim_ideal.expectation_value(observable)

    return np.array(expectations), ideal
