# Point Kinetics

Point kinetics is the simplest time-dependent model of nuclear systems. It treats the entire reactor or assembly as a single "point" with uniform properties, tracking only the total neutron population and energy release.

## The Point Kinetics Equations

The fundamental equation for neutron population is:

$$\frac{dN}{dt} = \alpha N$$

where α is the **Rossi alpha** (inverse reactor period):

$$\alpha = \frac{k_{eff} - 1}{\Lambda}$$

with:
- k_eff: effective multiplication factor
- Λ: prompt neutron generation time

For a supercritical system (k > 1), α > 0 and N grows exponentially.

## Including Delayed Neutrons

The full point kinetics equations with delayed neutrons are:

$$\frac{dN}{dt} = \frac{\rho - \beta}{\Lambda}N + \sum_i \lambda_i C_i$$

$$\frac{dC_i}{dt} = \frac{\beta_i}{\Lambda}N - \lambda_i C_i$$

where:
- ρ = (k-1)/k: reactivity
- β: total delayed neutron fraction
- β_i: delayed fraction for group i
- λ_i: decay constant for group i
- C_i: delayed neutron precursor concentration

For prompt supercritical systems (ρ > β), the delayed neutrons become irrelevant—the system responds on the prompt timescale.

## Rossi Alpha Calculation

The Python code computes α from diffusion theory:

$$\alpha = v D (B_m^2 - B_g^2)$$

where:
- v: average neutron velocity (~18 × 10⁶ m/s for fission spectrum)
- D: diffusion coefficient
- B_m²: material buckling
- B_g²: geometric buckling

```python
def inverse_bomb_period(data, radius, compression, extrapolation_factor=0.71045):
    params = compute_diffusion_params(data, compression)
    
    delta = extrapolation_factor * params["lambda_tr"]
    R_prime = radius + delta
    B_g_squared = (np.pi / R_prime) ** 2
    
    alpha = NEUTRON_VELOCITY * params["D"] * (params["B_m_squared"] - B_g_squared)
    
    return alpha, params
```

## Time Scales

For prompt supercritical Pu-239:

| Quantity | Value | Meaning |
|----------|-------|---------|
| Λ | ~10 ns | Prompt generation time |
| α | ~10⁸ /s | Inverse period |
| τ = 1/α | ~10 ns | e-folding time |

The neutron population doubles every τ×ln(2) ≈ 7 ns. In 1 μs (~100 doublings), the population increases by 10³⁰.

## Coupled Physics

Real explosions involve coupled physics:

1. **Neutronics**: dN/dt = αN
2. **Energy deposition**: dE/dt = Σ_f × N × v × ε_f
3. **Hydrodynamics**: Core expansion from pressure
4. **Feedback**: Expansion reduces density → reduces α

The Python code couples these:

```python
def point_kinetics_rhs(t, state, data, initial_radius, compression, core_mass):
    R, R_dot, N, E_fission, E_kinetic = state
    
    # Update compression from radius
    current_compression = (initial_radius / R) ** 3 * compression
    
    # Neutronics
    alpha, params = inverse_bomb_period(data, R, current_compression)
    dN_dt = alpha * N
    
    # Fission power
    power = params["Sigma_f"] * N * NEUTRON_VELOCITY * ENERGY_PER_FISSION_J
    dE_fission_dt = power
    
    # Hydrodynamics
    E_radiation = max(0.0, E_fission - E_kinetic)
    volume = (4/3) * np.pi * R**3
    pressure = E_radiation / (3 * volume)
    
    # Shell acceleration
    dR_dot_dt = 4 * np.pi * R**2 * pressure / m_eff
    dE_kinetic_dt = 4 * np.pi * R**2 * pressure * R_dot
    
    return [R_dot, dR_dot_dt, dN_dt, dE_fission_dt, dE_kinetic_dt]
```

## The Shutdown Process

A nuclear explosion terminates when expansion makes the system subcritical (α < 0).

The sequence:
1. **Initiation**: Neutron triggers chain reaction
2. **Exponential growth**: N doubles every ~7 ns
3. **Energy deposition**: Fission heats material
4. **Pressure buildup**: Radiation pressure dominates
5. **Expansion**: Core expands at km/s
6. **Shutdown**: Density drops → α becomes negative
7. **Disassembly**: System blows apart

The yield depends on how long the system stays supercritical during expansion.

## Typical Results

For 6.2 kg Pu-239 at 2.5× compression:

```
Initial conditions:
  Radius: 4.17 cm
  Compression: 2.50×
  Alpha: 2.3×10⁸ /s
  k_eff: 1.36

Results:
  Final yield: 12.4 kt
  Efficiency: 3.2%
  Peak temperature: 4.2×10⁸ K
  Peak pressure: 85 Gbar
  Time: 0.42 μs
```

## Limitations of Point Kinetics

Point kinetics assumes:
- Spatially uniform flux (poor near boundaries)
- Separable time dependence (flux shape doesn't change)
- Single effective parameters (no energy structure)

These assumptions break down for:
- Highly non-uniform systems
- Strongly localized criticality
- Multi-region configurations

For these cases, we need spatial kinetics (2D diffusion or Monte Carlo).

## Comparison with Monte Carlo

| Aspect | Point Kinetics | Monte Carlo |
|--------|---------------|-------------|
| Computes | α, k_eff, yield | k_eff, flux distribution |
| Speed | Very fast (~seconds) | Slower (~minutes) |
| Time-dependence | Natural | Requires special methods |
| Geometry | Sphere (analytical) | Any shape |
| Physics | One-group average | Continuous energy |

Point kinetics is excellent for:
- Parameter studies
- Understanding trends
- Coupled physics simulations
- Yield estimates

Monte Carlo is better for:
- Accurate k_eff
- Complex geometries
- Detailed physics
- Benchmark validation

## Code Structure

The Python point kinetics module:

```
point_kinetics.py
├── NuclearData          # Material properties
├── compute_macroscopic_cross_sections()
├── compute_diffusion_params()
├── compute_critical_radius()
├── compute_critical_mass()
├── inverse_bomb_period()   # α calculation
├── simulate_explosion()    # Time-dependent solver
└── plot_results()
```

Key parameters from ENDF/B-VIII.0:
- σ_f = 1.800 barns (fission)
- σ_s = 5.935 barns (total scattering)
- σ_c = 0.065 barns (capture)
- ν = 3.165 (neutrons/fission)

## Summary

Point kinetics provides:
- Fast parameter exploration
- Physical intuition
- Yield estimates
- Compression scaling studies

It complements Monte Carlo by enabling rapid design space exploration before committing to expensive detailed calculations.
