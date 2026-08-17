# Discrepancy Notes: Temme et al. 2017 vs Our Reproduction

This document investigates **where and why** our reproduction differs from the original paper results, and what those differences reveal about the method.

---

## Discrepancy 1: ZNE Overshoot at High Noise (High Confidence)

**Our Finding**: At base noise $\lambda \geq 0.05$ with three-level Richardson extrapolation (scales 1, 3, 5), ZNE **increases** the error rather than reducing it.

**Quantitative Result**:

| Base $\lambda$ | Raw Error $|\langle O\rangle_\lambda - \langle O\rangle_0|$ | ZNE Error | Better? |
| :---: | :---: | :---: | :---: |
| 0.01 | ~0.026 | ~0.003 | ✅ ZNE improves ~9× |
| 0.02 | ~0.051 | ~0.010 | ✅ ZNE improves ~5× |
| 0.05 | ~0.129 | ~0.238 | ❌ ZNE overshoots! |

**Physical Explanation**: Richardson extrapolation relies on the Taylor expansion $\langle O\rangle_\lambda \approx \langle O\rangle_0 + a_1\lambda + a_2\lambda^2 + ...$, which is only valid for small $\lambda$. When the maximum amplified noise reaches $\lambda_{\text{max}} = 5 \times 0.05 = 0.25$, higher-order terms dominate and the linear approximation fails catastrophically.

**What the paper does**: Temme et al. 2017 explicitly work in the regime $\lambda \ll 1$ and note that the method is designed for near-term NISQ hardware with "small" coherent error rates. Current IBM hardware has $\lambda \approx 0.001$–$0.01$ per 2-qubit gate, well within the working regime.

**Lesson**: ZNE is not a global error correction. It is a perturbative error *mitigation* valid only in the low-noise regime.

---

## Discrepancy 2: Absolute Scale of $\langle ZZ \rangle$ Decay (Minor)

**Our Finding**: Our decay rate at $\lambda = 0.05$ gives $\langle ZZ \rangle \approx 0.87$, consistent with the analytic prediction $(1 - 4\lambda/3)^d$ for $d=2$ noisy gates.

**Paper result**: Uses IBM hardware data (2017) showing similar but slightly different decay due to additional hardware noise sources (crosstalk, $T_1/T_2$ relaxation) not in our pure depolarizing model.

**Status**: Our reproduction matches the theoretical prediction. The deviation from hardware data is expected and acknowledged.

---

## Discrepancy 3: Exponential Extrapolation Performance

**Our Finding**: At $\lambda = 0.05$, exponential extrapolation performs *worse* than Richardson for depolarizing noise, with errors approximately 2× larger.

**Explanation**: Exponential extrapolation is physically motivated for amplitude-damping ($T_1$ relaxation) channels where $\langle O\rangle_\lambda \propto e^{-B\lambda}$. For depolarizing noise (which has a linear Taylor expansion to first order), polynomial/Richardson is a better fit.

**Recommendation**: Match extrapolation method to the dominant physical noise channel.
