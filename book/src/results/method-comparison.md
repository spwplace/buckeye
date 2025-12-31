# Method Comparison

This project implements the same nuclear physics problem using two fundamentally different approaches:

1. **Julia**: Monte Carlo neutron transport
2. **Python**: Deterministic diffusion theory (1-group, 2-group, and 2D spatial)

Let's compare them across multiple dimensions.

## Algorithmic Differences

### Monte Carlo (Julia)

```
For each generation:
    For each neutron:
        While alive:
            Sample collision distance (exponential)
            Sample reaction type (discrete)
            If fission: bank new neutrons
            If scatter: update energy and direction
            If capture or escape: kill neutron
    Estimate k from fission/source ratio
    Normalize fission bank for next generation
```

**Philosophy**: Simulate physics exactly, average over many trials.

### Diffusion Theory (Python)

```
Build diffusion coefficient D = 1/(3Σ_tr)
Compute material buckling B_m² = (νΣ_f - Σ_a)/D
Compute geometric buckling B_g² = (π/R')²
Criticality when B_m² = B_g²
```

**Philosophy**: Approximate physics, solve analytically or with deterministic numerics.

## Comparison Table

| Aspect | Monte Carlo (Julia) | Diffusion (Python) |
|--------|--------------------|--------------------|
| **Physics model** | Transport equation | P1 approximation |
| **Energy treatment** | Continuous | 1-2 groups |
| **Geometry** | Any shape | Simple shapes or mesh |
| **Angular dependence** | Explicit | Integrated out |
| **Boundary conditions** | Exact (vacuum) | Extrapolated |
| **Solution method** | Random sampling | Linear algebra |
| **Output** | k_eff ± σ | k_eff (no uncertainty) |
| **Convergence** | Statistical (1/√N) | Numerical (O(h²)) |

## Cross Section Treatment

### Monte Carlo
- Tabulated on energy grid (24 points from 10 eV to 20 MeV)
- Log-log interpolation
- Separate fission, elastic, inelastic, capture
- Watt spectrum for fission neutrons

### Diffusion
- Fission-spectrum averaged (single values)
- Transport-corrected: Σ_tr = Σ_t - μ̄Σ_s
- Grouped data (1 or 2 groups)

## Critical Mass Results

For bare δ-phase Pu-239:

| Method | Critical Radius | Critical Mass | Error vs Benchmark |
|--------|----------------|---------------|-------------------|
| Monte Carlo | 6.4 cm | 16.3 kg | +6% |
| Point kinetics | 4.9 cm | 10.0 kg | -39% |
| 2D FVM | 5.0 cm | 10.4 kg | -36% |
| Benchmark (Jezebel) | 6.385 cm | 16.2 kg | — |

**Key observation**: Monte Carlo is 6% high due to simplified cross sections. Diffusion methods are ~35-40% low due to the diffusion approximation breaking down for this optically thin system.

## Why Diffusion Fails Here

The diffusion approximation requires:
1. Many mean free paths across the system
2. Weak absorption
3. Nearly isotropic scattering

For bare Pu-239:
- System is ~1.5 mean free paths across (not many)
- Absorption is significant (Σ_a/Σ_t ≈ 0.25)
- Scattering is isotropic ✓

The approximation is at its limits, explaining the ~35% error in critical mass.

## Time-Dependent Comparison

### Point Kinetics (Python)

Solves coupled ODEs:
- dN/dt = αN
- dE/dt = power from N
- dR/dt = hydrodynamic expansion

Results for 6.2 kg Pu at 2.5× compression:
- Yield: ~12 kt
- Efficiency: ~3%
- Time: ~0.4 μs

### Monte Carlo Kinetics (Not implemented)

Would require:
- Time-dependent source iteration
- Delayed neutron precursors
- Much more computation

For prompt supercritical systems, point kinetics is adequate.

## Code Complexity

### Julia Monte Carlo (~1800 lines)
```
NuclearFootball.jl    68 lines   Module definition
materials.jl         137 lines   Nuclide/material data
cross_sections.jl    564 lines   XS lookup, interpolation, Watt spectrum
geometry.jl          482 lines   Sphere, ellipsoid, football shapes
particle.jl          339 lines   Transport physics
criticality.jl       385 lines   Power iteration, k-eigenvalue
analysis.jl          362 lines   Football analysis, validation
```

### Python Diffusion (~2500 lines)
```
nuclear_physics.py   844 lines   Full coupled simulation
point_kinetics.py    467 lines   Simpler standalone version
diffusion_2d.py      394 lines   2D axisymmetric FVM
spatial_neutronics.py 516 lines  More advanced 2D FVM
two_group_physics.py  606 lines  Two-group treatment
```

## Runtime Comparison

For k_eff calculation on the Jezebel benchmark:

| Method | Runtime | k_eff | Uncertainty |
|--------|---------|-------|-------------|
| Monte Carlo (10k particles, 150 gen) | ~120 s | 1.057 | ±0.003 |
| Point kinetics | <1 s | 0.987* | N/A |
| 2D FVM (60×120 mesh) | ~10 s | 1.004 | N/A |

*Point kinetics gives α and k_∞, not k_eff directly.

## When to Use Each

### Use Monte Carlo When:
- Accuracy matters (k_eff to ~1%)
- Geometry is complex
- Energy dependence is important
- Validating against benchmarks

### Use Diffusion When:
- Speed matters
- Exploring parameter space
- Need physical intuition
- Coupled physics (hydrodynamics)
- Analytical solutions available

## The Football Problem

### Monte Carlo Approach (Julia)
- Model exact superellipsoid shape
- Track neutrons through pointed ends
- Calculate k_eff = 1.21 ± 0.004
- Find critical mass ~17 kg in football shape

### Diffusion Approach (Python)
- Approximate as ellipsoid or sphere
- Use analytical buckling for prolate spheroid
- Underestimate critical mass significantly
- Still useful for compression/yield studies

## Complementary Strengths

The two approaches are **complementary**, not competitive:

1. **Python diffusion** for:
   - Fast parameter sweeps
   - Compression effects
   - Yield vs mass studies
   - Physical understanding

2. **Julia Monte Carlo** for:
   - Final k_eff calculation
   - Football geometry
   - Benchmark validation
   - Publication-quality results

## Summary

| Question | Best Tool | Why |
|----------|----------|-----|
| "What's k_eff for this football?" | Monte Carlo | Exact geometry |
| "How does yield scale with compression?" | Point kinetics | Fast + coupled physics |
| "What's the flux shape?" | 2D FVM | Spatial resolution |
| "Is this configuration critical?" | Either | Both give qualitative answer |
| "What's the critical mass to 5%?" | Monte Carlo | Accuracy |

Both codebases contribute to answering our central question: **Yes, a plutonium football would be supercritical**—whether calculated with Monte Carlo (k ≈ 1.21) or estimated with diffusion theory (deeply supercritical based on mass >> critical mass).
