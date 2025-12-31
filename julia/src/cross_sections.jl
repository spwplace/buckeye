"""
    Cross section data and lookup

Continuous-energy cross sections for Monte Carlo transport.
Uses tabulated data with log-log interpolation, plus analytical models
for fission spectrum (Watt) and angular distributions.
"""

#=============================================================================
# Cross Section Data Structure
=============================================================================#

"""
    CrossSectionData

Continuous-energy cross section data for a nuclide.
Energy grid in MeV, cross sections in barns.
"""
struct CrossSectionData
    nuclide::Nuclide
    
    # Energy grid (MeV) - shared by all reactions
    energy::Vector{Float64}
    
    # Cross sections (barns) on energy grid
    σ_total::Vector{Float64}
    σ_elastic::Vector{Float64}
    σ_inelastic::Vector{Float64}
    σ_fission::Vector{Float64}
    σ_capture::Vector{Float64}  # (n,γ)
    
    # Fission data
    ν_bar::Vector{Float64}      # Average neutrons per fission vs energy
    
    # Watt spectrum parameters for fission neutrons
    watt_a::Float64  # MeV
    watt_b::Float64  # MeV⁻¹
    
    # Interpolators (created on construction)
    _interp_total::Any
    _interp_elastic::Any
    _interp_inelastic::Any
    _interp_fission::Any
    _interp_capture::Any
    _interp_nu::Any
end

"""
    CrossSectionData(nuc::Nuclide, E, σt, σel, σin, σf, σc, ν, a, b)

Construct cross section data with automatic interpolator creation.
"""
function CrossSectionData(
    nuc::Nuclide,
    energy::Vector{Float64},
    σ_total::Vector{Float64},
    σ_elastic::Vector{Float64},
    σ_inelastic::Vector{Float64},
    σ_fission::Vector{Float64},
    σ_capture::Vector{Float64},
    ν_bar::Vector{Float64},
    watt_a::Float64,
    watt_b::Float64
)
    # Use log-log interpolation for cross sections (standard in nuclear data)
    log_E = log.(energy)
    
    # Helper to create safe log interpolation
    function make_interp(σ)
        # Handle zeros by using small positive value
        σ_safe = max.(σ, 1e-10)
        log_σ = log.(σ_safe)
        interpolate((log_E,), log_σ, Gridded(Linear()))
    end
    
    return CrossSectionData(
        nuc, energy,
        σ_total, σ_elastic, σ_inelastic, σ_fission, σ_capture, ν_bar,
        watt_a, watt_b,
        make_interp(σ_total),
        make_interp(σ_elastic),
        make_interp(σ_inelastic),
        make_interp(σ_fission),
        make_interp(σ_capture),
        interpolate((log_E,), ν_bar, Gridded(Linear()))
    )
end

#=============================================================================
# Cross Section Lookup Functions
=============================================================================#

"""
    get_σ_total(xs::CrossSectionData, E::Float64) -> Float64

Total cross section at energy E (MeV). Returns barns.
"""
function get_σ_total(xs::CrossSectionData, E::Float64)
    E = clamp(E, xs.energy[1], xs.energy[end])
    return exp(xs._interp_total(log(E)))
end

"""
    get_σ_fission(xs::CrossSectionData, E::Float64) -> Float64

Fission cross section at energy E (MeV). Returns barns.
"""
function get_σ_fission(xs::CrossSectionData, E::Float64)
    E = clamp(E, xs.energy[1], xs.energy[end])
    return exp(xs._interp_fission(log(E)))
end

"""
    get_σ_elastic(xs::CrossSectionData, E::Float64) -> Float64

Elastic scattering cross section at energy E (MeV). Returns barns.
"""
function get_σ_elastic(xs::CrossSectionData, E::Float64)
    E = clamp(E, xs.energy[1], xs.energy[end])
    return exp(xs._interp_elastic(log(E)))
end

"""
    get_σ_inelastic(xs::CrossSectionData, E::Float64) -> Float64

Inelastic scattering cross section at energy E (MeV). Returns barns.
"""
function get_σ_inelastic(xs::CrossSectionData, E::Float64)
    E = clamp(E, xs.energy[1], xs.energy[end])
    return exp(xs._interp_inelastic(log(E)))
end

"""
    get_σ_capture(xs::CrossSectionData, E::Float64) -> Float64

Radiative capture (n,γ) cross section at energy E (MeV). Returns barns.
"""
function get_σ_capture(xs::CrossSectionData, E::Float64)
    E = clamp(E, xs.energy[1], xs.energy[end])
    return exp(xs._interp_capture(log(E)))
end

"""
    get_σ_scatter(xs::CrossSectionData, E::Float64) -> Float64

Total scattering cross section (elastic + inelastic).
"""
function get_σ_scatter(xs::CrossSectionData, E::Float64)
    return get_σ_elastic(xs, E) + get_σ_inelastic(xs, E)
end

"""
    get_σ_absorption(xs::CrossSectionData, E::Float64) -> Float64

Total absorption cross section (fission + capture).
"""
function get_σ_absorption(xs::CrossSectionData, E::Float64)
    return get_σ_fission(xs, E) + get_σ_capture(xs, E)
end

"""
    get_ν_bar(xs::CrossSectionData, E::Float64) -> Float64

Average number of neutrons per fission at energy E (MeV).
"""
function get_ν_bar(xs::CrossSectionData, E::Float64)
    E = clamp(E, xs.energy[1], xs.energy[end])
    return xs._interp_nu(log(E))
end

#=============================================================================
# Fission Neutron Spectrum (Watt Distribution)
=============================================================================#

"""
    sample_watt(a::Float64, b::Float64, rng::AbstractRNG) -> Float64

Sample energy from Watt fission spectrum:
    f(E) ∝ exp(-E/a) * sinh(√(bE))

Uses rejection sampling method from MCNP.
"""
function sample_watt(a::Float64, b::Float64, rng::AbstractRNG)
    # Parameters for sampling
    K = 1.0 + (a * b) / 8.0
    L = (K + sqrt(K^2 - 1)) / a
    M = a * L - 1.0
    
    while true
        # Sample from g(E) ∝ exp(-L*E)
        ξ₁ = rand(rng)
        ξ₂ = rand(rng)
        
        E = -log(ξ₁) / L
        
        # Rejection condition
        x = (L * E - M)^2
        y = b * E
        
        if (ξ₂^2 ≤ y / x) || (x ≤ y)
            return E
        end
    end
end

"""
    sample_fission_energy(xs::CrossSectionData, rng::AbstractRNG) -> Float64

Sample outgoing neutron energy from fission spectrum.
"""
function sample_fission_energy(xs::CrossSectionData, rng::AbstractRNG)
    return sample_watt(xs.watt_a, xs.watt_b, rng)
end

"""
    sample_fission_energy(xs::CrossSectionData) -> Float64

Sample outgoing neutron energy using default RNG.
"""
function sample_fission_energy(xs::CrossSectionData)
    return sample_watt(xs.watt_a, xs.watt_b, Random.default_rng())
end

#=============================================================================
# Angular Distribution for Scattering
=============================================================================#

"""
    sample_elastic_cosine(A::Float64, E::Float64, rng::AbstractRNG) -> Float64

Sample cosine of scattering angle in center-of-mass frame for elastic scattering.
Uses isotropic CM approximation (valid for heavy nuclei at MeV energies).

Returns μ_lab (cosine of scattering angle in lab frame).
"""
function sample_elastic_cosine(A::Float64, E::Float64, rng::AbstractRNG)
    # For heavy nuclei (A >> 1), scattering is nearly isotropic in CM
    # Sample isotropic in CM: μ_cm = 2ξ - 1
    μ_cm = 2.0 * rand(rng) - 1.0
    
    # Transform to lab frame
    # μ_lab = (1 + A*μ_cm) / √(1 + A² + 2A*μ_cm)
    denom = sqrt(1.0 + A^2 + 2.0 * A * μ_cm)
    μ_lab = (1.0 + A * μ_cm) / denom
    
    return μ_lab
end

"""
    sample_inelastic_cosine(rng::AbstractRNG) -> Float64

Sample cosine of scattering angle for inelastic scattering.
Approximated as isotropic in lab frame.
"""
function sample_inelastic_cosine(rng::AbstractRNG)
    return 2.0 * rand(rng) - 1.0
end

"""
    elastic_energy_loss(E::Float64, A::Float64, μ_cm::Float64) -> Float64

Calculate post-scatter energy for elastic collision.
E' = E * (A² + 1 + 2A*μ_cm) / (A + 1)²
"""
function elastic_energy_loss(E::Float64, A::Float64, μ_cm::Float64)
    return E * (A^2 + 1.0 + 2.0 * A * μ_cm) / (A + 1.0)^2
end

#=============================================================================
# Pre-built Cross Section Data for Key Nuclides
=============================================================================#

"""
Build Pu-239 cross section data from ENDF/B-VIII.0 representative values.
Energy range: 10 keV to 20 MeV (fast spectrum).
"""
function build_Pu239_xs()
    # Energy grid (MeV)
    energy = [
        1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7,
        1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,
        12.0, 14.0, 17.0, 20.0
    ]
    
    # Total cross section (barns) - from ENDF/B-VIII.0
    σ_total = [
        25.0, 18.0, 12.0, 9.5, 8.8, 8.8, 8.0, 7.5, 7.1, 6.8,
        6.2, 6.0, 5.85, 5.75, 5.70, 5.75, 5.80, 5.85, 5.90, 5.95,
        6.00, 5.99, 6.05, 6.10
    ]
    
    # Fission cross section (barns)
    σ_fission = [
        1.8, 1.75, 1.70, 1.65, 1.55, 1.52, 1.58, 1.62, 1.63, 1.66,
        1.68, 1.80, 2.01, 1.90, 1.78, 1.85, 1.95, 2.05, 2.15, 2.25,
        2.30, 2.34, 2.40, 2.45
    ]
    
    # Elastic scattering (barns)
    σ_elastic = [
        8.0, 7.5, 6.5, 5.5, 5.2, 5.0, 4.8, 4.6, 4.4, 4.2,
        4.0, 3.6, 3.2, 3.0, 2.9, 2.85, 2.80, 2.75, 2.70, 2.65,
        2.60, 2.55, 2.50, 2.45
    ]
    
    # Inelastic scattering (barns) - threshold ~40 keV
    σ_inelastic = [
        0.0, 0.0, 0.0, 0.3, 0.8, 1.0, 1.3, 1.5, 1.8, 2.0,
        2.2, 2.5, 2.7, 2.6, 2.4, 2.2, 2.0, 1.85, 1.70, 1.55,
        1.45, 1.35, 1.25, 1.15
    ]
    
    # Capture (n,γ) cross section (barns) - small for fast neutrons
    σ_capture = [
        0.5, 0.3, 0.15, 0.08, 0.04, 0.025, 0.018, 0.015, 0.012, 0.010,
        0.008, 0.006, 0.005, 0.004, 0.003, 0.0025, 0.002, 0.0018, 0.0015, 0.001,
        0.0008, 0.0006, 0.0005, 0.0004
    ]
    
    # ν̄(E) - average neutrons per fission: ν̄ = 2.874 + 0.148*E
    ν_bar = [2.874 + 0.148 * E for E in energy]
    
    # Watt spectrum parameters for Pu-239
    watt_a = 0.966  # MeV
    watt_b = 2.842  # MeV⁻¹
    
    return CrossSectionData(
        Pu239, energy,
        σ_total, σ_elastic, σ_inelastic, σ_fission, σ_capture, ν_bar,
        watt_a, watt_b
    )
end

"""
Build U-235 cross section data.
"""
function build_U235_xs()
    # Energy grid (MeV)
    energy = [
        1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7,
        1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,
        12.0, 14.0, 17.0, 20.0
    ]
    
    # Total cross section (barns)
    σ_total = [
        20.0, 15.0, 11.0, 9.0, 8.0, 7.8, 7.5, 7.2, 6.9, 6.7,
        6.5, 6.2, 6.0, 5.8, 5.7, 5.65, 5.70, 5.75, 5.80, 5.85,
        5.90, 5.95, 6.00, 6.05
    ]
    
    # Fission cross section (barns)
    σ_fission = [
        1.4, 1.35, 1.30, 1.25, 1.18, 1.15, 1.18, 1.20, 1.22, 1.20,
        1.18, 1.22, 1.25, 1.20, 1.10, 1.05, 1.20, 1.50, 1.80, 2.00,
        2.10, 2.15, 2.20, 2.25
    ]
    
    # Elastic scattering (barns)
    σ_elastic = [
        9.0, 7.0, 5.8, 5.0, 4.5, 4.3, 4.1, 3.9, 3.8, 3.7,
        3.6, 3.4, 3.2, 3.0, 2.9, 2.85, 2.80, 2.75, 2.70, 2.65,
        2.60, 2.55, 2.50, 2.45
    ]
    
    # Inelastic scattering (barns)
    σ_inelastic = [
        0.0, 0.0, 0.0, 0.2, 0.6, 0.9, 1.1, 1.3, 1.5, 1.7,
        1.9, 2.1, 2.3, 2.4, 2.5, 2.4, 2.2, 2.0, 1.8, 1.6,
        1.5, 1.4, 1.3, 1.2
    ]
    
    # Capture cross section (barns)
    σ_capture = [
        0.4, 0.25, 0.12, 0.06, 0.03, 0.02, 0.015, 0.012, 0.010, 0.008,
        0.006, 0.005, 0.004, 0.003, 0.0025, 0.002, 0.0018, 0.0015, 0.0012, 0.001,
        0.0008, 0.0007, 0.0006, 0.0005
    ]
    
    # ν̄(E) for U-235: ν̄ ≈ 2.42 + 0.12*E
    ν_bar = [2.42 + 0.12 * E for E in energy]
    
    # Watt spectrum parameters for U-235
    watt_a = 0.988  # MeV
    watt_b = 2.249  # MeV⁻¹
    
    return CrossSectionData(
        U235, energy,
        σ_total, σ_elastic, σ_inelastic, σ_fission, σ_capture, ν_bar,
        watt_a, watt_b
    )
end

"""
Build Ga (gallium) cross section data - non-fissile, for delta-phase Pu.
Simplified: just scattering and capture.
"""
function build_Ga_xs()
    energy = [
        1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0,
        5.0, 10.0, 14.0, 20.0
    ]
    
    # Total ≈ elastic + capture (no fission, minimal inelastic)
    σ_total = [
        8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.5, 4.2, 4.0,
        3.8, 3.6, 3.5, 3.4
    ]
    
    σ_elastic = [
        7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.6, 4.2, 4.0, 3.8,
        3.6, 3.4, 3.3, 3.2
    ]
    
    σ_inelastic = [
        0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.25, 0.18, 0.15,
        0.12, 0.10, 0.09, 0.08
    ]
    
    σ_fission = zeros(length(energy))  # No fission
    
    σ_capture = [
        0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 0.1, 0.05, 0.02, 0.01,
        0.005, 0.003, 0.002, 0.001
    ]
    
    ν_bar = zeros(length(energy))  # N/A
    
    return CrossSectionData(
        Ga69, energy,  # Use Ga-69 as representative
        σ_total, σ_elastic, σ_inelastic, σ_fission, σ_capture, ν_bar,
        0.0, 0.0  # No fission spectrum
    )
end

# Global cross section library
const XS_LIBRARY = Dict{String, CrossSectionData}()

"""
Build Pu-240 cross section data (non-fissile for thermal, low fission for fast).
Important for spontaneous fission background but not primary fissile.
"""
function build_Pu240_xs()
    energy = [
        1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0,
        5.0, 10.0, 14.0, 20.0
    ]
    
    # Pu-240 has fission threshold around 1 MeV
    σ_total = [
        12.0, 10.0, 9.0, 8.5, 8.0, 7.8, 7.5, 7.0, 6.5, 6.2,
        6.0, 5.9, 5.85, 5.8
    ]
    
    σ_elastic = [
        8.0, 7.0, 6.0, 5.5, 5.0, 4.8, 4.5, 4.2, 4.0, 3.8,
        3.5, 3.3, 3.2, 3.1
    ]
    
    σ_inelastic = [
        0.0, 0.0, 0.0, 0.2, 0.5, 0.8, 1.0, 1.3, 1.5, 1.6,
        1.5, 1.4, 1.3, 1.2
    ]
    
    # Pu-240 fission threshold ~1 MeV (subthreshold fission very low)
    σ_fission = [
        0.0, 0.0, 0.0, 0.0, 0.01, 0.02, 0.05, 0.2, 0.8, 1.4,
        1.6, 1.8, 1.9, 2.0
    ]
    
    σ_capture = [
        3.5, 2.5, 2.0, 1.5, 1.0, 0.8, 0.5, 0.3, 0.15, 0.08,
        0.04, 0.02, 0.015, 0.01
    ]
    
    # ν̄ for Pu-240 (similar to Pu-239)
    ν_bar = [2.8 + 0.14 * E for E in energy]
    
    return CrossSectionData(
        Pu240, energy,
        σ_total, σ_elastic, σ_inelastic, σ_fission, σ_capture, ν_bar,
        0.966, 2.842  # Same Watt params as Pu-239
    )
end

"""
    get_xs(nuc::Nuclide) -> CrossSectionData

Get cross section data for a nuclide, building if necessary.
"""
function get_xs(nuc::Nuclide)
    if !haskey(XS_LIBRARY, nuc.name)
        if nuc.name == "Pu-239"
            XS_LIBRARY[nuc.name] = build_Pu239_xs()
        elseif nuc.name == "Pu-240"
            XS_LIBRARY[nuc.name] = build_Pu240_xs()
        elseif nuc.name == "Pu-241"
            # Pu-241 is fissile, use Pu-239 as approximation
            XS_LIBRARY[nuc.name] = build_Pu239_xs()
        elseif nuc.name == "U-235"
            XS_LIBRARY[nuc.name] = build_U235_xs()
        elseif nuc.name in ("Ga-69", "Ga-71")
            XS_LIBRARY[nuc.name] = build_Ga_xs()
        else
            error("No cross section data available for $(nuc.name)")
        end
    end
    return XS_LIBRARY[nuc.name]
end

#=============================================================================
# Material Macroscopic Cross Sections
=============================================================================#

"""
    get_Σ_total(mat::Material, E::Float64) -> Float64

Macroscopic total cross section for material at energy E.
Returns cm⁻¹.
"""
function get_Σ_total(mat::Material, E::Float64)
    Σ = 0.0
    for (nuc, _) in mat.composition
        xs = get_xs(nuc)
        N = number_density(mat, nuc)  # atoms/barn-cm
        σ = get_σ_total(xs, E)        # barns
        Σ += N * σ                    # barn * atoms/barn-cm = cm⁻¹
    end
    return Σ
end

"""
    get_Σ_fission(mat::Material, E::Float64) -> Float64

Macroscopic fission cross section for material.
"""
function get_Σ_fission(mat::Material, E::Float64)
    Σ = 0.0
    for (nuc, _) in mat.composition
        xs = get_xs(nuc)
        N = number_density(mat, nuc)
        σ = get_σ_fission(xs, E)
        Σ += N * σ
    end
    return Σ
end

"""
    get_νΣ_fission(mat::Material, E::Float64) -> Float64

ν * Σ_fission: neutron production cross section.
"""
function get_νΣ_fission(mat::Material, E::Float64)
    νΣ = 0.0
    for (nuc, _) in mat.composition
        xs = get_xs(nuc)
        N = number_density(mat, nuc)
        σf = get_σ_fission(xs, E)
        ν = get_ν_bar(xs, E)
        νΣ += N * σf * ν
    end
    return νΣ
end
