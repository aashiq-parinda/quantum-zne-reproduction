"""Temme et al. 2017 Figure Reproduction Analysis.

Reproduces the key quantitative results from:
    Temme, K., Bravyi, S., Gambetta, J. M. (2017).
    Error mitigation for short-depth quantum circuits.
    Physical Review Letters, 119(18), 180509.
    arXiv: 1612.02058

Figures We Reproduce:
    Fig 1: ⟨ZZ⟩ expectation vs noise strength λ (linear decay)
    Fig 2: ZNE vs raw error rate as a function of circuit depth
    Fig 3: Mitigation performance across different extrapolation methods

Test circuit: Bell state preparation circuit on 2 qubits
    |ψ⟩ = H ⊗ I → CNOT → |Φ+⟩ = (|00⟩ + |11⟩)/√2
    Observable: Z⊗Z, ideal value = ⟨ZZ⟩_0 = 1.0
"""

import numpy as np
from qzne.simulation.noisy_simulator import (
    DensityMatrixSimulator, build_zz_observable, H, X, Z
)
from qzne.simulation.zne_pipeline import zne_workflow, scan_noise_levels


def bell_state_circuit(sim: DensityMatrixSimulator, n_folds: int = 0) -> None:
    """Bell state preparation: H(0) → CNOT(0,1).

    With gate folding: each gate applied (2*n_folds+1) times to amplify noise.
    """
    # Apply H to qubit 0 (with folding)
    for _ in range(2 * n_folds + 1):
        sim.apply_unitary(H, 0)
        sim.apply_noise(0)
    # Apply CNOT(0→1) (with folding)
    for _ in range(2 * n_folds + 1):
        sim.apply_cnot(0, 1)
        sim.apply_noise(0)
        sim.apply_noise(1)


def ghz_circuit(sim: DensityMatrixSimulator, n_folds: int = 0) -> None:
    """GHZ state: H(0) → CNOT(0,1) → CNOT(1,2) on 3 qubits."""
    for _ in range(2 * n_folds + 1):
        sim.apply_unitary(H, 0)
        sim.apply_noise(0)
    for _ in range(2 * n_folds + 1):
        sim.apply_cnot(0, 1)
        sim.apply_noise(0)
        sim.apply_noise(1)
    for _ in range(2 * n_folds + 1):
        sim.apply_cnot(1, 2)
        sim.apply_noise(1)
        sim.apply_noise(2)


def reproduce_figure1() -> dict:
    """Reproduce Fig 1: ⟨ZZ⟩ vs noise strength λ for Bell state.

    Expectation: Linear decay from ⟨ZZ⟩=1.0 at λ=0 to ⟨ZZ⟩=0 at λ=3/4.
    Paper shows the observable decays as ⟨ZZ⟩_λ ≈ (1 - 4λ/3)^d
    where d = circuit depth (number of noisy gates).
    """
    lambdas = np.linspace(0.0, 0.20, 30)
    obs = build_zz_observable(2, 0, 1)
    expectations, ideal = scan_noise_levels(bell_state_circuit, 2, obs, lambdas)

    return {
        "lambdas": lambdas,
        "expectations_zz": expectations,
        "ideal_value": ideal,
        "description": "Fig 1: <ZZ> vs depolarizing noise rate λ (Bell state, 2 qubits)",
    }


def reproduce_figure2() -> dict:
    """Reproduce Fig 2: ZNE error vs raw error across noise levels.

    Compares:
    - Raw noisy ⟨O⟩_λ  (no mitigation)
    - ZNE estimate from Richardson extrapolation at {λ, 3λ, 5λ}
    """
    results = []
    base_lambdas = np.linspace(0.01, 0.12, 15)
    obs = build_zz_observable(2, 0, 1)

    for base_lam in base_lambdas:
        res = zne_workflow(
            circuit_fn=bell_state_circuit,
            n_qubits=2,
            base_lambda=base_lam,
            noise_scales=[1, 3, 5],
            observable=obs,
            method="richardson",
        )
        results.append({
            "base_lambda": base_lam,
            "raw_error": res["raw_error"],
            "zne_error": res["zne_error"],
            "mitigation_factor": res["mitigation_factor"],
            "zne_estimate": res["zne_estimate"],
            "ideal": res["ideal_value"],
        })

    return {
        "data": results,
        "description": "Fig 2: Raw error vs ZNE error at varying base noise rates",
    }


def reproduce_figure3() -> dict:
    """Reproduce Fig 3: Comparison of ZNE extrapolation methods.

    Compares: linear, Richardson (3-point), polynomial (degree 2), exponential
    at noise rate λ = 0.05.
    """
    base_lam = 0.05
    obs = build_zz_observable(2, 0, 1)
    methods = ["linear", "richardson", "polynomial", "exponential"]
    method_results = {}

    for method in methods:
        scales = [1, 3] if method == "linear" else [1, 3, 5]
        res = zne_workflow(
            circuit_fn=bell_state_circuit,
            n_qubits=2,
            base_lambda=base_lam,
            noise_scales=scales,
            observable=obs,
            method=method,
        )
        method_results[method] = {
            "zne_estimate": res["zne_estimate"],
            "zne_error": res["zne_error"],
            "mitigation_factor": res["mitigation_factor"],
            "raw_error": res["raw_error"],
            "ideal": res["ideal_value"],
        }

    return {
        "base_lambda": base_lam,
        "methods": method_results,
        "description": "Fig 3: ZNE extrapolation method comparison at λ=0.05",
    }


if __name__ == "__main__":
    print("=" * 70)
    print(" 📊 TEMME et al. 2017 — ZNE FIGURE REPRODUCTION")
    print("=" * 70)

    fig1 = reproduce_figure1()
    print(f"\n  Fig 1 — Noise Scan Summary:")
    print(f"  Ideal ⟨ZZ⟩_0 = {fig1['ideal_value']:.6f}")
    print(f"  At λ=0.05: ⟨ZZ⟩ = {fig1['expectations_zz'][7]:.6f}")
    print(f"  At λ=0.10: ⟨ZZ⟩ = {fig1['expectations_zz'][14]:.6f}")

    fig2 = reproduce_figure2()
    data = fig2["data"]
    avg_mitigation = np.mean([d["mitigation_factor"] for d in data])
    print(f"\n  Fig 2 — ZNE Mitigation Performance:")
    print(f"  Average mitigation factor: {avg_mitigation:.2f}×")
    print(f"  Best case reduction: {max(d['mitigation_factor'] for d in data):.2f}×")

    fig3 = reproduce_figure3()
    print(f"\n  Fig 3 — Method Comparison (λ=0.05):")
    for method, d in fig3["methods"].items():
        print(f"    {method:12s}: ZNE={d['zne_estimate']:+.6f}  error={d['zne_error']:.6f}  factor={d['mitigation_factor']:.2f}×")
