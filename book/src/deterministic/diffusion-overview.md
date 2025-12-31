# Diffusion Theory Overview

While the Julia code uses Monte Carlo neutron transport, the Python codebase takes a fundamentally different approach: **deterministic diffusion theory**. This chapter introduces diffusion theory and explains when it's appropriate versus Monte Carlo.

## The Diffusion Approximation

Instead of tracking individual neutrons, diffusion theory treats neutrons as a continuous field that "diffuses" through matter like heat or chemical concentration.

### From Transport to Diffusion

The full neutron transport equation is an integro-differential equation in 7 dimensions (3 space + 2 angle + 1 energy + 1 time). The diffusion approximation simplifies this by assuming:

1. **Near-isotropic angular flux**: The neutron angular distribution is only weakly anisotropic
2. **Gradual spatial variation**: Flux changes slowly compared to mean free path
3. **No strong absorbers**: Absorption is weak relative to scattering

Under these assumptions, the angular dependence can be integrated out, yielding Fick's Law:

$$\mathbf{J} = -D \nabla \phi$$

where:
- **J** is the neutron current (net flow vector)
- **D** is the diffusion coefficient
- **φ** is the scalar neutron flux

### The Diffusion Equation

Combining Fick's Law with neutron balance gives the diffusion equation:

$$\frac{1}{v}\frac{\partial \phi}{\partial t} = D\nabla^2\phi - \Sigma_a\phi + \nu\Sigma_f\phi + S$$

In steady state without external sources:

$$D\nabla^2\phi + (\nu\Sigma_f - \Sigma_a)\phi = 0$$

This is an eigenvalue problem—we seek the flux distribution φ that satisfies this equation.

## The Diffusion Coefficient

The diffusion coefficient is:

$$D = \frac{1}{3\Sigma_{tr}}$$

where Σ_tr is the **transport cross section**:

$$\Sigma_{tr} = \Sigma_t - \bar{\mu}\Sigma_s$$

The term \\(\bar{\mu} = \langle \cos\theta \rangle\\) accounts for preferential forward scattering. For heavy nuclei:

$$\bar{\mu} \approx \frac{2}{3A}$$

For Pu-239 (A=239): \\(\bar{\mu} \approx 0.0028\\), so Σ_tr ≈ Σ_t.

## Material Buckling

The key parameter in diffusion theory is the **material buckling**:

$$B_m^2 = \frac{\nu\Sigma_f - \Sigma_a}{D}$$

This represents the "curvature" the flux would have in an infinite medium. Positive B² means the material is intrinsically supercritical.

For k_∞ > 1:
$$B_m^2 = \frac{k_\infty - 1}{L^2}$$

where L² = D/Σ_a is the diffusion length squared.

## Geometric Buckling

In a finite geometry, the flux must satisfy boundary conditions. For a sphere of radius R with zero flux at the extrapolated boundary R' = R + δ:

$$\phi(r) = \frac{A}{r}\sin(Br)$$

where B must satisfy:

$$BR' = \pi \quad \Rightarrow \quad B = \frac{\pi}{R'}$$

The **geometric buckling** is:

$$B_g^2 = \left(\frac{\pi}{R'}\right)^2$$

## The Criticality Condition

A system is critical when geometric buckling equals material buckling:

$$B_g^2 = B_m^2$$

This gives the critical radius:

$$R_c = \frac{\pi}{\sqrt{B_m^2}} - \delta$$

where δ is the extrapolation distance (typically 0.71λ_tr for vacuum boundaries).

## Extrapolation Distance

At a vacuum boundary, the flux doesn't go to zero exactly at the surface—it extrapolates to zero a short distance beyond. The Milne problem gives:

$$\delta = 0.7104 \times \lambda_{tr}$$

This correction is crucial for accurate critical mass predictions. The Python code implements both:
- Standard Milne value (0.7104)
- Trinity-calibrated value (0.633) to match experimental yields

## When Diffusion Theory Works

Diffusion theory is accurate when:
- System is many mean free paths across
- Absorption is moderate (c = Σ_s/Σ_t > 0.5)
- No strong spatial gradients
- Far from boundaries

For a bare Pu-239 sphere:
- λ_tr ≈ 4 cm
- Critical radius ≈ 6.4 cm
- System is ~1.5 mean free paths across

This is marginal for diffusion theory! The approximation works reasonably well but has ~10-15% error on critical mass compared to transport calculations.

## Diffusion Theory in the Python Code

The Python implementation (`nuclear_physics.py`) uses:

```python
def compute_diffusion_params(data: NuclearData, compression: float) -> DiffusionParameters:
    xs = compute_macroscopic_xs(data, compression)
    
    lambda_tr = 1.0 / xs.Sigma_tr
    D = lambda_tr / 3.0
    
    # Buckling-dependent correction for absorbing media
    correction = 2.0 * D * xs.Sigma_a / (5.0 + 3.0 * D * xs.Sigma_a)
    D_tr = D * (1.0 - correction)
    
    B_m_sq = (data.nu_prompt * xs.Sigma_f - xs.Sigma_a) / D_tr
    L_sq = D_tr / xs.Sigma_a
    k_inf = data.nu_prompt * xs.Sigma_f / xs.Sigma_a
    
    return DiffusionParameters(D=D, D_tr=D_tr, B_m_sq=B_m_sq, ...)
```

The code includes transport corrections that improve accuracy for absorbing media.

## Diffusion vs Monte Carlo

| Aspect | Diffusion Theory | Monte Carlo |
|--------|-----------------|-------------|
| **Basis** | Continuum approximation | Individual particle tracking |
| **Computation** | Solve differential equations | Random sampling |
| **Speed** | Fast (seconds) | Slow (minutes to hours) |
| **Accuracy** | ~10-15% error on critical mass | ~1% with good XS data |
| **Geometry** | Analytical for simple shapes | Any geometry |
| **Energy** | One or few groups | Continuous |
| **Near boundaries** | Poor (approximation breaks down) | Exact |
| **Time-dependence** | Natural extension | More complex |

## Summary

Diffusion theory provides:
- Analytical solutions for simple geometries
- Physical insight into criticality
- Fast parameter studies
- Reasonable (~15%) accuracy for bulk properties

But it struggles with:
- Strong absorbers
- Thin systems (few mean free paths)
- Complex geometries
- Accurate boundary treatment

Our Python code uses diffusion theory for:
- Critical mass estimation
- Time-dependent kinetics (point kinetics)
- Parameter studies (compression effects)

The Julia code uses Monte Carlo for:
- Accurate k_eff calculation
- Complex geometry (football shape)
- Benchmark validation
- Continuous energy physics

Both approaches give valuable insights; they complement rather than replace each other.
