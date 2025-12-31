# Criticality and the Multiplication Factor

The central question in nuclear criticality is simple: does the neutron population grow, shrink, or stay constant? The answer is captured in a single number called **k**, the effective multiplication factor.

## The Multiplication Factor

Consider a generation of N neutrons in a fissile assembly. After they undergo various reactions (fission, scattering, absorption, leakage), a new generation of neutrons is produced. The ratio of successive generations defines k:

$$k = \frac{N_{n+1}}{N_n} = \frac{\text{neutrons in generation } n+1}{\text{neutrons in generation } n}$$

Three regimes exist:

| Condition | Name | Behavior |
|-----------|------|----------|
| k < 1 | Subcritical | Neutron population decreases exponentially |
| k = 1 | Critical | Neutron population remains constant |
| k > 1 | Supercritical | Neutron population increases exponentially |

## The Six-Factor Formula

For thermal reactors, k can be decomposed into the **six-factor formula**:

$$k = \eta f p \epsilon P_{NL}^T P_{NL}^F$$

where:

- **η (eta)**: Neutrons produced per thermal neutron absorbed in fuel (~2.1 for Pu-239)
- **f**: Thermal utilization—fraction of absorptions that occur in fuel
- **p**: Resonance escape probability—fraction of neutrons that avoid resonance capture while slowing down
- **ε**: Fast fission factor—boost from fast fissions in U-238
- **P_NL^T**: Thermal non-leakage probability
- **P_NL^F**: Fast non-leakage probability

For our bare plutonium football, this simplifies considerably:
- Pure fuel: f = 1
- No moderator: p = 1, ε ≈ 1
- No thermal neutrons: η replaced by fast-spectrum average

The remaining physics is dominated by the non-leakage probability.

## k-effective vs k-infinite

Two versions of k are commonly used:

### k_∞ (k-infinite)
The multiplication factor for an infinitely large system with no leakage:

$$k_\infty = \frac{\text{neutron production}}{\text{absorption}}$$

For Pu-239, using energy-averaged cross sections:
$$k_\infty = \frac{\bar{\nu} \Sigma_f}{\Sigma_a} = \frac{\bar{\nu} \sigma_f}{\sigma_f + \sigma_c}$$

With ν̄ ≈ 2.9, σ_f ≈ 1.7 b, and σ_c ≈ 0.01 b:
$$k_\infty \approx \frac{2.9 \times 1.7}{1.7 + 0.01} \approx 2.88$$

This is greater than 2—pure Pu-239 has ample criticality margin.

### k_eff (k-effective)
The actual multiplication factor including leakage:

$$k_{eff} = k_\infty P_{NL}$$

where P_NL is the non-leakage probability. For a finite assembly:

$$P_{NL} = \frac{\text{neutrons absorbed}}{\text{neutrons absorbed + escaped}}$$

The non-leakage probability depends on geometry, size, and material. A sphere minimizes surface area (and thus leakage) for a given volume.

## Criticality Calculation

For a given geometry and material, finding k_eff is an **eigenvalue problem**. The neutron transport equation can be written:

$$\hat{L}\phi = \frac{1}{k}\hat{F}\phi$$

where:
- $\hat{L}$ is the loss operator (absorption + leakage)
- $\hat{F}$ is the fission production operator
- $\phi$ is the neutron flux
- $k$ is the eigenvalue we seek

The physical interpretation: the fission source must be scaled by 1/k to maintain steady state.

### Power Iteration

The standard numerical method is **power iteration**:

1. Guess an initial fission source distribution S₀(r)
2. Solve for the flux: $\hat{L}\phi_n = S_n$
3. Calculate the new source: $S_{n+1} = \hat{F}\phi_n$
4. Estimate k: $k_n = \frac{\int S_{n+1}}{\int S_n}$
5. Normalize: $S_{n+1} \leftarrow S_{n+1} / k_n$
6. Repeat until converged

Monte Carlo implements this by tracking neutrons generation by generation.

## Reactivity

The **reactivity** ρ measures the departure from critical:

$$\rho = \frac{k - 1}{k}$$

| k | ρ | Status |
|---|---|--------|
| 0.95 | -5.26% | Subcritical |
| 0.99 | -1.01% | Slightly subcritical |
| 1.00 | 0 | Critical |
| 1.01 | +0.99% | Slightly supercritical |
| 1.05 | +4.76% | Supercritical |

Reactivity is often expressed in "dollars" where 1 dollar = β (the delayed neutron fraction). For Pu-239:
- 1 dollar = 0.22% = 220 pcm
- Prompt critical occurs at ρ = β (1 dollar)

Our football calculation will show ρ ≈ +17%, or about 75 dollars supercritical—deeply in the prompt supercritical regime.

## Prompt vs Delayed Criticality

### Delayed Critical (0 < ρ < β)
The system is critical including delayed neutrons. The neutron population grows slowly, with a period of seconds to minutes. Reactors operate in this regime.

### Prompt Critical (ρ = β)
The system is critical on prompt neutrons alone. At this point, the delayed neutrons become irrelevant for population growth. The reactor period drops to milliseconds.

### Prompt Supercritical (ρ > β)
Each prompt neutron generation produces more than one prompt neutron. Exponential growth occurs on a timescale of:

$$T = \frac{l}{k_p - 1}$$

where l is the prompt neutron generation time (~10 ns for fast systems) and k_p is the prompt multiplication factor.

For our football with k ≈ 1.2:
$$T \approx \frac{10^{-8}}{0.2} = 50 \text{ ns}$$

The neutron population doubles every 50 nanoseconds—this is a nuclear explosion.

## Critical Mass

The **critical mass** is the minimum mass of fissile material needed to achieve k = 1 for a given geometry.

Critical mass depends strongly on:

1. **Material**: Pu-239 < U-233 < U-235
2. **Density**: Higher density → smaller critical mass
3. **Shape**: Sphere is optimal
4. **Reflector**: Surrounding material can reduce critical mass by 30-50%
5. **Enrichment**: Higher fissile fraction → lower critical mass

### Bare Critical Masses

For bare (unreflected) spheres:

| Material | Density (g/cm³) | Critical Mass (kg) | Critical Radius (cm) |
|----------|-----------------|-------------------|---------------------|
| Pu-239 (α-phase) | 19.86 | 10.0 | 4.9 |
| Pu-239 (δ-phase) | 15.8 | 16.3 | 6.4 |
| U-233 | 18.8 | 16.4 | 5.8 |
| U-235 | 18.7 | 52.0 | 8.7 |

Delta-phase plutonium has 63% higher critical mass than alpha-phase due to its lower density. But delta-phase is used in weapons because alpha-phase is brittle and pyrophoric.

### The Football

Our NFL football holds about 22 kg of delta-phase plutonium—about 35% above critical mass. Even accounting for the suboptimal shape (higher leakage than a sphere), this is supercritical.

## Geometry Effects

The non-leakage probability depends on the **buckling**, which quantifies how curved the flux distribution is:

$$B^2 = \frac{\nabla^2 \phi}{\phi}$$

For a sphere of radius R:
$$B^2 = \left(\frac{\pi}{R}\right)^2$$

For an ellipsoid with semi-axes a, b, c:
$$B^2 = \frac{\pi^2}{a^2} + \frac{\pi^2}{b^2} + \frac{\pi^2}{c^2}$$

The non-leakage probability is approximately:
$$P_{NL} \approx \frac{1}{1 + L^2 B^2}$$

where L is the neutron diffusion length.

A sphere minimizes buckling (and thus maximizes P_NL) for a given volume. Any other shape has higher leakage and requires more mass to go critical.

For our football shape:
- Volume equivalent to a sphere of radius ~8.2 cm
- But actual half-length is 14 cm with pointed ends
- Effective buckling is ~15% higher than equivalent sphere
- Critical mass is ~15% higher than spherical configuration

## Summary

| Parameter | Value for Pu-239 Football |
|-----------|---------------------------|
| k_∞ | ~2.9 |
| P_NL | ~0.42 |
| k_eff | ~1.2 |
| ρ | ~17% |
| Status | Prompt supercritical |

The football is deeply supercritical. In the next chapter, we'll see how Monte Carlo methods calculate k_eff by simulating individual neutrons.
