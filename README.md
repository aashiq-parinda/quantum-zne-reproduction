# Reproduction of Temme et al. 2017: Zero Noise Extrapolation (ZNE)

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Tests: 21/21 Passed](https://img.shields.io/badge/tests-21%2F21%20passing-brightgreen)](#-testing)
[![DOI: 10.5281/zenodo.21979332](https://zenodo.org/badge/DOI/10.5281/zenodo.21979332.svg)](https://doi.org/10.5281/zenodo.21979332)

This repository reproduces the foundational quantum error mitigation paper:
> **Temme, K., Bravyi, S., & Gambetta, J. M. (2017).** *Error mitigation for short-depth quantum circuits.* Physical Review Letters, 119(18), 180509. [arXiv:1612.02058](https://arxiv.org/abs/1612.02058) | [Zenodo Reproduction Record](https://doi.org/10.5281/zenodo.21979332).

It is built from mathematical first principles using a custom density-matrix quantum simulator with depolarizing noise channels, gate folding, and Richardson extrapolation.

---

## 📐 Mathematical Core

### 1. Physical Model & Noise Amplification (Gate Folding)
Quantum circuits subject to depolarizing noise undergo channel transformations:
$$\mathcal{E}_\lambda(\rho) = (1-\lambda)\rho + \frac{\lambda}{3}(X\rho X + Y\rho Y + Z\rho Z)$$

Gate folding scales noise by replacing unitaries $U \rightarrow U(U^\dagger U)^n$, scaling noise by odd integer factors $c = 2n + 1 \in \{1, 3, 5, \dots\}$ without altering the ideal logic.

### 2. Richardson Extrapolation
Expectation values $\langle O \rangle_\lambda$ at noise rates $\lambda_1, \dots, \lambda_k$ are extrapolated to zero noise ($\lambda \to 0$):
$$\langle O \rangle_0 \approx \sum_{k=1}^n c_k \langle O \rangle_{\lambda_k} \quad \text{where } c_k = \prod_{j \neq k} \frac{\lambda_j}{\lambda_j - \lambda_k}$$

---

## 📊 Reproduction Results & Figures

| Figure | Paper Target | Reproduction Finding | Status |
| :--- | :--- | :--- | :---: |
| **Fig 1** | $\langle ZZ \rangle$ vs noise rate $\lambda$ | Confirmed linear decay matching theoretical $(1 - 4\lambda/3)^d$ | ✅ Reproduced |
| **Fig 2** | Raw error vs ZNE error | Confirmed ~5–10× error reduction in low-noise regime ($\lambda \le 0.02$) | ✅ Reproduced |
| **Fig 3** | Method comparison | Richardson outperforms linear 2-point extrapolation | ✅ Reproduced |

---

## ⚠️ Documented Discrepancy & Regimes

As documented in [`docs/DISCREPANCY_NOTES.md`](docs/DISCREPANCY_NOTES.md), Richardson extrapolation **overshoots** at higher noise levels ($\lambda \ge 0.05$). Because 3-scale folding amplifies effective noise to $5\lambda = 0.25$, the higher-order terms dominate, violating the small-$\lambda$ Taylor expansion assumption. This confirms that ZNE is valid specifically in the NISQ low-noise regime ($\lambda \sim 0.001 - 0.01$).

---

## 💻 Quickstart & Verification

```bash
# Clone & set up
git clone https://github.com/aashiq-parinda/quantum-zne-reproduction.git
cd quantum-zne-reproduction
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# Run full reproduction demonstration
python example.py

# Run test suite (21 unit tests)
pytest tests/ -v
```

---

## 📜 Structure & References

- `src/qzne/core/`: Noise channels, gate folding, and Richardson extrapolation math.
- `src/qzne/simulation/`: Density matrix simulator and pipeline.
- `src/qzne/analysis/`: Reproduction code for Figures 1, 2, and 3.
- `docs/PAPER_SUMMARY.md`: Summary of paper methodology and assumptions.
- `docs/DISCREPANCY_NOTES.md`: Investigation of extrapolation limits.
