# Cross Section Physics

Cross sections quantify the probability of nuclear reactions. Getting them right is essential for accurate Monte Carlo simulation.

## Microscopic Cross Sections

The **microscopic cross section** σ(E) is the effective target area presented by a nucleus to a neutron of energy E. It has units of area (barns, where 1 b = 10⁻²⁴ cm²).

Conceptually:
$$\sigma = \frac{\text{reaction probability per target nucleus}}{\text{neutron fluence}}$$

A nucleus is not a hard sphere—cross sections can be much larger (resonance capture) or smaller (high-energy scattering) than the geometric area.

## Types of Cross Sections

For neutron reactions, we need several partial cross sections:

| Symbol | Reaction | Meaning |
|--------|----------|---------|
| σ_t | Total | Any interaction |
| σ_el | Elastic scatter | (n,n) |
| σ_in | Inelastic scatter | (n,n') |
| σ_f | Fission | (n,f) |
| σ_γ | Radiative capture | (n,γ) |
| σ_a | Absorption | σ_f + σ_γ |

The total is the sum of all partials:
$$\sigma_t = \sigma_{el} + \sigma_{in} + \sigma_f + \sigma_\gamma + ...$$

## Energy Dependence

Cross sections vary dramatically with energy. Three regimes:

### Thermal Region (E < 1 eV)

Many cross sections follow the **1/v law**:
$$\sigma(E) \propto \frac{1}{v} \propto \frac{1}{\sqrt{E}}$$

Physically: slower neutrons spend more time near the nucleus.

### Resonance Region (1 eV - 10 keV)

Sharp peaks occur when the neutron energy matches excited states of the compound nucleus. The **Breit-Wigner formula**:

$$\sigma(E) = \sigma_0 \frac{\Gamma^2/4}{(E-E_0)^2 + \Gamma^2/4}$$

where:
- E₀ is the resonance energy
- Γ is the total width (inverse lifetime)
- σ₀ is the peak cross section

Heavy nuclei like Pu-239 have thousands of resonances.

### Fast Region (> 10 keV)

Resonances overlap and average out. Cross sections become relatively smooth:
- Fission: ~1-2 barns for Pu-239
- Elastic scatter: ~3-5 barns
- Capture: Falls to ~0.01 barns

## Pu-239 Cross Section Data

Our code uses tabulated values derived from ENDF/B-VIII.0:

| E (MeV) | σ_t (b) | σ_f (b) | σ_el (b) | σ_in (b) | σ_γ (b) |
|---------|---------|---------|----------|----------|---------|
| 0.01 | 9.5 | 1.65 | 5.5 | 0.3 | 0.08 |
| 0.1 | 8.8 | 1.52 | 5.0 | 1.0 | 0.025 |
| 1.0 | 6.2 | 1.68 | 4.0 | 2.2 | 0.008 |
| 5.0 | 5.75 | 1.85 | 2.85 | 2.2 | 0.0025 |
| 10.0 | 5.95 | 2.25 | 2.65 | 1.55 | 0.001 |

Key observations:
- Fission cross section relatively flat in fast region (~1.5-2.5 b)
- Elastic scatter dominates at low energy
- Inelastic scatter has threshold (~40 keV)
- Capture is negligible for fast neutrons

## Interpolation

Between tabulated points, we interpolate. Standard practice is **log-log interpolation**:

$$\log(\sigma) = \log(\sigma_1) + \frac{\log(E) - \log(E_1)}{\log(E_2) - \log(E_1)}[\log(\sigma_2) - \log(\sigma_1)]$$

This works well because cross sections often follow power laws.

```julia
function get_σ(xs::CrossSectionData, E::Float64)
    # Clamp to data range
    E = clamp(E, xs.energy[1], xs.energy[end])
    # Log-log interpolation
    return exp(xs._interp(log(E)))
end
```

## Macroscopic Cross Sections

For materials with number density N (atoms/cm³), the **macroscopic cross section** is:
$$\Sigma = N \sigma$$

Units: cm⁻¹ (probability of interaction per cm traveled).

For a material with multiple nuclides:
$$\Sigma_t = \sum_i N_i \sigma_{t,i}$$

### Number Density Calculation

For material with density ρ and atomic mass A:
$$N = \frac{\rho N_A}{A}$$

For Pu-239 at 15.8 g/cm³:
$$N = \frac{15.8 \times 6.022 \times 10^{23}}{239} = 3.98 \times 10^{22} \text{ atoms/cm}^3$$

Converting to atoms/barn-cm (conventional units):
$$N = 0.0398 \text{ atoms/barn-cm}$$

## ν̄: Neutrons per Fission

The average number of neutrons produced per fission varies with incident energy:
$$\bar{\nu}(E) = \bar{\nu}_0 + \alpha E$$

For Pu-239:
- ν̄₀ = 2.874
- α = 0.148 MeV⁻¹

At E = 1 MeV: ν̄ ≈ 3.02

The quantity νΣf (neutron production cross section) is crucial for criticality:
$$\nu\Sigma_f = N \bar{\nu} \sigma_f$$

## The Fission Spectrum

Fission neutrons are born with a characteristic energy distribution. The **Watt spectrum** is:

$$\chi(E) = C e^{-E/a} \sinh(\sqrt{bE})$$

where:
- a = 0.966 MeV for Pu-239
- b = 2.842 MeV⁻¹ for Pu-239
- C is a normalization constant

Properties:
- Peak energy: ~0.7 MeV
- Average energy: ~2.0 MeV
- Extends from 0 to ~15 MeV

### Sampling the Watt Spectrum

We use rejection sampling (the MCNP algorithm):

```julia
function sample_watt(a, b, rng)
    K = 1.0 + (a * b) / 8.0
    L = (K + sqrt(K^2 - 1)) / a
    M = a * L - 1.0
    
    while true
        ξ₁ = rand(rng)
        ξ₂ = rand(rng)
        
        E = -log(ξ₁) / L
        
        x = (L * E - M)^2
        y = b * E
        
        if ξ₂^2 ≤ y / x || x ≤ y
            return E
        end
    end
end
```

This samples from an exponential envelope and accepts with probability proportional to sinh(√bE).

## Scattering Angular Distributions

### Elastic Scattering

For heavy nuclei, elastic scattering in the center-of-mass frame is approximately isotropic:
$$P(\mu_{cm}) = \frac{1}{2}, \quad \mu_{cm} \in [-1, 1]$$

Transformation to lab frame:
$$\mu_{lab} = \frac{1 + A\mu_{cm}}{\sqrt{A^2 + 2A\mu_{cm} + 1}}$$

For A = 239 (Pu), even μ_cm = -1 gives μ_lab ≈ 0.004 (essentially isotropic in lab too).

### Energy Loss in Elastic Scatter

The post-scatter energy:
$$E' = E \cdot \frac{A^2 + 2A\mu_{cm} + 1}{(A + 1)^2}$$

Define:
$$\alpha = \left(\frac{A-1}{A+1}\right)^2$$

For Pu-239: α = 0.983. A neutron can lose at most 1.7% of its energy per collision.

The average logarithmic energy loss:
$$\xi = 1 + \frac{\alpha \ln \alpha}{1 - \alpha} \approx \frac{2}{A + 2/3}$$

For Pu-239: ξ ≈ 0.0084.

## Cross Section Library Structure

Our code organizes cross section data as:

```julia
struct CrossSectionData
    nuclide::Nuclide
    
    # Energy grid (MeV)
    energy::Vector{Float64}
    
    # Cross sections (barns)
    σ_total::Vector{Float64}
    σ_elastic::Vector{Float64}
    σ_inelastic::Vector{Float64}
    σ_fission::Vector{Float64}
    σ_capture::Vector{Float64}
    
    # Fission data
    ν_bar::Vector{Float64}
    watt_a::Float64
    watt_b::Float64
    
    # Interpolators
    _interp_total::Interpolation
    ...
end
```

Interpolators are built on construction for fast lookup during transport.

## ENDF Data Format

Production codes read from ENDF (Evaluated Nuclear Data File):
- MF=3: Cross sections
- MF=4: Angular distributions  
- MF=5: Energy distributions
- MF=6: Product energy-angle distributions

Our simplified tables capture the essential physics without the full complexity.

## Summary

| Quantity | Pu-239 Value | Units |
|----------|--------------|-------|
| σ_f (1 MeV) | 1.68 | barns |
| σ_t (1 MeV) | 6.2 | barns |
| ν̄ (1 MeV) | 3.02 | neutrons/fission |
| χ average | 2.0 | MeV |
| α (elastic) | 0.983 | dimensionless |
| λ (mfp) | 4.0 | cm |

These numbers fully characterize fast neutron behavior in plutonium.
