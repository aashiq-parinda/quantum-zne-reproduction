"""Density Matrix Noisy Quantum Circuit Simulator.

Executes quantum circuits under realistic noise models by propagating
the full density matrix ρ rather than a statevector.

For N qubits: ρ ∈ ℂ^{2^N × 2^N}, initialized as |0...0⟩⟨0...0|.

Gate application: ρ → UρU†
After-gate noise: ρ → ε_λ(ρ) = Σ_k K_k ρ K_k†

This simulator enables us to:
1. Compute exact noisy expectation values ⟨O⟩_λ = Tr(O ε_λ(ρ))
2. Scan ⟨O⟩_λ vs λ to reproduce Fig 1 of Temme et al. 2017
3. Validate ZNE extrapolation against ground truth ⟨O⟩_0
"""

import numpy as np
from typing import List, Callable, Dict, Any
from qzne.core.noise_models import apply_kraus_channel, depolarizing_kraus


# Standard Pauli gates
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)


class DensityMatrixSimulator:
    """N-qubit density matrix quantum circuit simulator with per-gate noise.

    Parameters
    ----------
    n_qubits : int
    noise_lambda : float
        Depolarizing noise strength per gate (0 = ideal).
    noise_type : str
        'depolarizing' (default), 'dephasing', 'none'
    """

    def __init__(
        self,
        n_qubits: int,
        noise_lambda: float = 0.0,
        noise_type: str = "depolarizing",
    ) -> None:
        self.n = n_qubits
        self.dim = 2 ** n_qubits
        self.noise_lambda = noise_lambda
        self.noise_type = noise_type
        self.reset()

    def reset(self) -> None:
        """Reset to |0...0⟩⟨0...0|."""
        self.rho = np.zeros((self.dim, self.dim), dtype=complex)
        self.rho[0, 0] = 1.0

    def _single_qubit_superop(
        self,
        gate: np.ndarray,
        qubit: int,
    ) -> np.ndarray:
        """Build 2^N × 2^N superoperator I⊗...⊗gate⊗...⊗I."""
        ops = [I2] * self.n
        ops[qubit] = gate
        M = ops[0]
        for op in ops[1:]:
            M = np.kron(M, op)
        return M

    def apply_unitary(self, gate: np.ndarray, qubit: int) -> None:
        """Apply single-qubit unitary U to qubit k: ρ → UρU†."""
        U = self._single_qubit_superop(gate, qubit)
        self.rho = U @ self.rho @ U.conj().T

    def apply_cnot(self, control: int, target: int) -> None:
        """Apply CNOT gate on 2-qubit subspace (control, target)."""
        if self.n < 2:
            raise ValueError("CNOT requires N ≥ 2 qubits")
        dim = self.dim
        P = np.zeros((dim, dim), dtype=complex)
        for idx in range(dim):
            bits = [(idx >> (self.n - 1 - b)) & 1 for b in range(self.n)]
            if bits[control] == 1:
                bits[target] ^= 1
            new_idx = sum(bits[b] << (self.n - 1 - b) for b in range(self.n))
            P[new_idx, idx] = 1
        self.rho = P @ self.rho @ P.conj().T

    def apply_noise(self, qubit: int) -> None:
        """Apply per-gate depolarizing noise to a single qubit.

        ε_λ(ρ) via tensor product embedding of the 2×2 channel into 2^N space.
        """
        if self.noise_lambda == 0.0 or self.noise_type == "none":
            return

        if self.noise_type == "depolarizing":
            kraus_ops_1q = depolarizing_kraus(self.noise_lambda)
        else:
            return

        # Build 2^N Kraus operators by tensor embedding
        rho_out = np.zeros_like(self.rho)
        ops = [I2] * self.n
        for K in kraus_ops_1q:
            ops[qubit] = K
            K_full = ops[0]
            for op in ops[1:]:
                K_full = np.kron(K_full, op)
            rho_out += K_full @ self.rho @ K_full.conj().T
        self.rho = rho_out

    def expectation_value(self, observable: np.ndarray) -> float:
        """Compute ⟨O⟩ = Tr(O ρ).

        Parameters
        ----------
        observable : np.ndarray, shape (2^N, 2^N) — Hermitian observable

        Returns
        -------
        float — real part of Tr(O ρ)
        """
        return float(np.real(np.trace(observable @ self.rho)))

    def purity(self) -> float:
        """Compute purity Tr(ρ²) ∈ [1/2^N, 1]."""
        return float(np.real(np.trace(self.rho @ self.rho)))


def build_zz_observable(n_qubits: int, q1: int = 0, q2: int = 1) -> np.ndarray:
    """Build Z⊗Z observable on qubits (q1, q2) in an N-qubit system."""
    ops = [I2] * n_qubits
    ops[q1] = Z
    ops[q2] = Z
    M = ops[0]
    for op in ops[1:]:
        M = np.kron(M, op)
    return M
