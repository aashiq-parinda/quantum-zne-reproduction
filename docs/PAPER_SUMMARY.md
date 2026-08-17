# Paper Summary: Temme et al. 2017 — Zero Noise Extrapolation

**Full Citation**: Temme, K., Bravyi, S., & Gambetta, J. M. (2017). *Error mitigation for short-depth quantum circuits*. Physical Review Letters, **119**(18), 180509. DOI: 10.1103/PhysRevLett.119.180509. arXiv: 1612.02058.

---

## Problem Statement

Near-term quantum devices ("NISQ") have significant gate errors and decoherence. Without error correction, every gate introduces noise:

$$\mathcal{E}_\lambda(\rho) = (1-\lambda)\rho + \frac{\lambda}{3}(X\rho X + Y\rho Y + Z\rho Z)$$

This corrupts expectation values of observables $O$:

$$\langle O \rangle_\lambda = \langle O \rangle_0 + \lambda \epsilon_1 + \lambda^2 \epsilon_2 + O(\lambda^3)$$

**Goal**: Estimate the noise-free value $\langle O \rangle_0$ without additional qubit resources.

---

## ZNE Method — Richardson Extrapolation

### Step 1: Noise Amplification via Gate Folding

Replace each gate $U$ with a folded sequence that has identical ideal action but amplified noise:

$$U \rightarrow U (U^\dagger U)^n \qquad \text{noise factor } c = 2n+1$$

This gives effective noise rates $\{c_1 \lambda, c_2 \lambda, c_3 \lambda, ...\}$ for scale factors $\{c_1, c_2, c_3, ...\} = \{1, 3, 5, ...\}$.

### Step 2: Measure Expectation Values

Execute circuit at each noise level and record $\langle O \rangle_{c_k \lambda}$.

### Step 3: Richardson Extrapolation

Compute the zero-noise estimate using Richardson's deferred approach to the limit:

$$\langle O \rangle_0 \approx \sum_{k=1}^{n} c_k \langle O \rangle_{\lambda_k}$$

where the Richardson coefficients are:

$$c_k = \prod_{j \neq k} \frac{\lambda_j}{\lambda_j - \lambda_k}$$

**Key property**: $\sum_k c_k = 1$, and this cancels all noise terms up to order $\lambda^{n-1}$.

---

## Key Assumptions

1. **Small noise** ($\lambda \ll 1$): Taylor expansion is valid. Method degrades at large $\lambda$.
2. **Markovian noise**: Noise is uncorrelated between gates.
3. **Gate-independent noise**: Same $\lambda$ applies to every gate.
4. **Locality**: Single-qubit depolarizing after each gate.

---

## What We Reproduce

| Figure | Original Result | Our Result | Status |
| :--- | :--- | :--- | :--- |
| Fig 1: $\langle ZZ \rangle$ vs $\lambda$ | Linear decay $\approx (1-4\lambda/3)^d$ | ✅ Linear decay confirmed | Reproduced |
| Fig 2: ZNE vs raw error | ZNE reduces error ~5-10× at low $\lambda$ | ✅ Confirmed at $\lambda \le 0.02$ | Reproduced |
| Fig 3: Method comparison | Richardson outperforms linear 2-point | ✅ Richardson has lower error | Reproduced |

See `DISCREPANCY_NOTES.md` for documented differences.
