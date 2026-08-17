# Paper Summary: Temme et al. 2017 — Zero Noise Extrapolation

**Full Citation**: Temme, K., Bravyi, S., & Gambetta, J. M. (2017). *Error mitigation for short-depth quantum circuits*. Physical Review Letters, **119**(18), 180509. DOI: 10.1103/PhysRevLett.119.180509. arXiv: 1612.02058.

---

## Problem Statement

Near-term quantum devices ("NISQ") have significant gate errors and decoherence. Without error correction, every gate introduces noise:

```
E_λ(ρ) = (1 - λ)ρ + (λ/3) · [XρX + YρY + ZρZ]
```

This corrupts expectation values of observables O:

```
⟨O⟩_λ = ⟨O⟩₀ + λ ε₁ + λ² ε₂ + O(λ³)
```

**Goal**: Estimate the noise-free value ⟨O⟩₀ without additional qubit resources.

---

## ZNE Method — Richardson Extrapolation

### Step 1: Noise Amplification via Gate Folding

Replace each gate U with a folded sequence that has identical ideal action but amplified noise:

```
U  →  U (U† U)^n     (noise scaling factor c = 2n + 1)
```

This gives effective noise rates {c₁λ, c₂λ, c₃λ, ...} for scale factors {c₁, c₂, c₃, ...} = {1, 3, 5, ...}.

### Step 2: Measure Expectation Values

Execute circuit at each noise level and record ⟨O⟩_{c_k λ}.

### Step 3: Richardson Extrapolation

Compute the zero-noise estimate using Richardson's deferred approach to the limit:

```
⟨O⟩₀ ≈ ∑_{k=1}^n c_k · ⟨O⟩_{λ_k}
```

where the Richardson coefficients are:

```
c_k = ∏_{j ≠ k} [ λ_j / (λ_j - λ_k) ]
```

**Key property**: ∑_k c_k = 1, and this cancels all noise terms up to order λ^(n-1).

---

## Key Assumptions

1. **Small noise** (λ ≪ 1): Taylor expansion is valid. Method degrades at large λ.
2. **Markovian noise**: Noise is uncorrelated between gates.
3. **Gate-independent noise**: Same λ applies to every gate.
4. **Locality**: Single-qubit depolarizing after each gate.

---

## What We Reproduce

| Figure | Original Result | Our Result | Status |
| :--- | :--- | :--- | :--- |
| Fig 1: ⟨ZZ⟩ vs λ | Linear decay ≈ (1 - 4λ/3)^d | ✅ Linear decay confirmed | Reproduced |
| Fig 2: ZNE vs raw error | ZNE reduces error ~5-10× at low λ | ✅ Confirmed at λ ≤ 0.02 | Reproduced |
| Fig 3: Method comparison | Richardson outperforms linear 2-point | ✅ Richardson has lower error | Reproduced |

See `DISCREPANCY_NOTES.md` for documented differences.
