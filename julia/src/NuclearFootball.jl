"""
    NuclearFootball

Monte Carlo neutron transport code to answer the critical question:
Can a nuclear football achieve criticality?

Uses continuous-energy Monte Carlo with k-eigenvalue power iteration.
"""
module NuclearFootball

using LinearAlgebra
using StaticArrays
using Random
using Statistics
using Printf
using Interpolations
using ProgressMeter

# Type aliases for clarity
const Vec3 = SVector{3, Float64}
const BARN_TO_CM2 = 1e-24  # 1 barn = 1e-24 cm²

# Physical constants
const NEUTRON_MASS_AMU = 1.008665
const AMU_TO_KG = 1.66054e-27
const MEV_TO_JOULE = 1.602176634e-13
const AVOGADRO = 6.02214076e23

# Include submodules in dependency order
include("materials.jl")
include("cross_sections.jl")
include("geometry.jl")
include("particle.jl")
include("criticality.jl")
include("analysis.jl")

# Main exports
export Vec3

# Materials
export Nuclide, Material
export Pu239, Pu240, U235, U238
export δ_phase_plutonium, α_phase_plutonium, weapons_grade_pu

# Geometry
export Geometry, Sphere, Ellipsoid, Football
export point_inside, distance_to_boundary, surface_normal
export volume, bounding_box

# Cross sections
export CrossSectionData, get_σ_total, get_σ_fission, get_σ_elastic
export get_σ_inelastic, get_σ_capture, get_ν_bar
export sample_fission_energy, sample_scatter_angle

# Particle transport
export Particle, SourceSite
export transport!, sample_collision_distance

# Criticality
export CriticalityResult, KEffEstimates
export run_criticality, find_critical_mass
export shannon_entropy

# Analysis
export analyze_football, compare_geometries, print_results, validate_jezebel

end # module
