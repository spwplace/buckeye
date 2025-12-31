"""
    Particle transport

Defines neutron particle state and transport mechanics.
"""

#=============================================================================
# Particle State
=============================================================================#

"""
    Particle

A neutron being transported through the geometry.
"""
mutable struct Particle
    r::Vec3           # Position (cm)
    u::Vec3           # Direction (unit vector)
    E::Float64        # Energy (MeV)
    weight::Float64   # Statistical weight
    alive::Bool       # Still being tracked?
end

"""
    Particle(r::Vec3, E::Float64)

Create a particle at position r with energy E, isotropic direction.
"""
function Particle(r::Vec3, E::Float64, rng::AbstractRNG)
    # Sample isotropic direction
    μ = 2*rand(rng) - 1  # cos(θ) uniform in [-1, 1]
    φ = 2π * rand(rng)
    
    sin_θ = sqrt(1 - μ^2)
    u = Vec3(sin_θ * cos(φ), sin_θ * sin(φ), μ)
    
    return Particle(r, u, E, 1.0, true)
end

"""
    SourceSite

A stored fission site for the fission bank.
"""
struct SourceSite
    r::Vec3           # Position
    E::Float64        # Energy (sampled from fission spectrum)
    weight::Float64   # Weight
end

#=============================================================================
# Direction Sampling
=============================================================================#

"""
    sample_isotropic_direction(rng::AbstractRNG) -> Vec3

Sample a uniformly distributed direction on the unit sphere.
"""
function sample_isotropic_direction(rng::AbstractRNG)
    μ = 2*rand(rng) - 1
    φ = 2π * rand(rng)
    sin_θ = sqrt(1 - μ^2)
    return Vec3(sin_θ * cos(φ), sin_θ * sin(φ), μ)
end

"""
    rotate_direction(u::Vec3, μ::Float64, rng::AbstractRNG) -> Vec3

Rotate direction u by scattering angle with cos(θ) = μ.
Azimuthal angle is uniformly random.
"""
function rotate_direction(u::Vec3, μ::Float64, rng::AbstractRNG)
    # Sample azimuthal angle
    φ = 2π * rand(rng)
    
    sin_θ = sqrt(max(0.0, 1 - μ^2))
    cos_φ = cos(φ)
    sin_φ = sin(φ)
    
    # If u is nearly aligned with z-axis, handle specially
    if abs(u[3]) > 0.9999
        u_new = Vec3(sin_θ * cos_φ, sin_θ * sin_φ, sign(u[3]) * μ)
    else
        # General rotation
        a = sqrt(1 - u[3]^2)
        u_new = Vec3(
            u[1]*μ + (sin_θ/a)*(u[1]*u[3]*cos_φ - u[2]*sin_φ),
            u[2]*μ + (sin_θ/a)*(u[2]*u[3]*cos_φ + u[1]*sin_φ),
            u[3]*μ - a*sin_θ*cos_φ
        )
    end
    
    # Renormalize to ensure unit vector
    return normalize(u_new)
end

#=============================================================================
# Collision Distance Sampling
=============================================================================#

"""
    sample_collision_distance(Σ_total::Float64, rng::AbstractRNG) -> Float64

Sample distance to next collision from exponential distribution.
d ~ -ln(ξ) / Σ_total
"""
function sample_collision_distance(Σ_total::Float64, rng::AbstractRNG)
    ξ = rand(rng)
    # Avoid log(0)
    while ξ == 0.0
        ξ = rand(rng)
    end
    return -log(ξ) / Σ_total
end

#=============================================================================
# Collision Physics
=============================================================================#

"""
    ReactionType

Types of neutron reactions.
"""
@enum ReactionType begin
    ELASTIC
    INELASTIC
    FISSION
    CAPTURE
    ESCAPE
end

"""
    sample_reaction(mat::Material, E::Float64, rng::AbstractRNG) -> ReactionType

Sample which reaction occurs at collision.
"""
function sample_reaction(mat::Material, E::Float64, rng::AbstractRNG)
    Σ_total = get_Σ_total(mat, E)
    
    # Accumulate partial cross sections
    Σ_elastic = 0.0
    Σ_inelastic = 0.0
    Σ_fission = 0.0
    Σ_capture = 0.0
    
    for (nuc, _) in mat.composition
        xs = get_xs(nuc)
        N = number_density(mat, nuc)
        Σ_elastic   += N * get_σ_elastic(xs, E)
        Σ_inelastic += N * get_σ_inelastic(xs, E)
        Σ_fission   += N * get_σ_fission(xs, E)
        Σ_capture   += N * get_σ_capture(xs, E)
    end
    
    # Sample reaction type
    ξ = rand(rng) * Σ_total
    
    if ξ < Σ_elastic
        return ELASTIC
    elseif ξ < Σ_elastic + Σ_inelastic
        return INELASTIC
    elseif ξ < Σ_elastic + Σ_inelastic + Σ_fission
        return FISSION
    else
        return CAPTURE
    end
end

"""
    get_dominant_fissile(mat::Material) -> Nuclide

Get the dominant fissile nuclide in a material (for fission spectrum sampling).
"""
function get_dominant_fissile(mat::Material)
    # Return the nuclide with highest fission cross section
    best_nuc = mat.composition[1][1]
    best_σf = 0.0
    
    for (nuc, _) in mat.composition
        xs = get_xs(nuc)
        σf = get_σ_fission(xs, 1.0)  # Check at 1 MeV
        if σf > best_σf
            best_σf = σf
            best_nuc = nuc
        end
    end
    
    return best_nuc
end

#=============================================================================
# Transport Step
=============================================================================#

"""
    transport_step!(p::Particle, geom::Geometry, mat::Material, 
                   fission_bank::Vector{SourceSite}, rng::AbstractRNG,
                   k_eff::Float64) -> Tuple{Float64, Float64}

Perform one transport step for a particle.
Returns (track_length_tally, collision_tally) contributions to k_eff.
"""
function transport_step!(
    p::Particle, 
    geom::Geometry, 
    mat::Material,
    fission_bank::Vector{SourceSite}, 
    rng::AbstractRNG,
    k_eff::Float64
)
    track_length_keff = 0.0
    collision_keff = 0.0
    
    # Get macroscopic cross section
    Σ_total = get_Σ_total(mat, p.E)
    
    # Sample collision distance
    d_collision = sample_collision_distance(Σ_total, rng)
    
    # Distance to boundary
    d_boundary = distance_to_boundary(geom, p.r, p.u)
    
    if d_collision < d_boundary
        # Collision occurs inside geometry
        
        # Track-length estimator contribution
        νΣf = get_νΣ_fission(mat, p.E)
        track_length_keff = p.weight * d_collision * νΣf
        
        # Move to collision site
        p.r = p.r + d_collision * p.u
        
        # Sample reaction
        reaction = sample_reaction(mat, p.E, rng)
        
        # Collision estimator contribution
        collision_keff = p.weight * νΣf / Σ_total
        
        if reaction == FISSION
            # Bank fission neutrons for next generation
            # Expected number: ν̄ * weight / k_eff
            nuc = get_dominant_fissile(mat)
            xs = get_xs(nuc)
            ν = get_ν_bar(xs, p.E)
            
            # Number to bank (stochastic rounding)
            expected = p.weight * ν / k_eff
            n_bank = floor(Int, expected)
            if rand(rng) < (expected - n_bank)
                n_bank += 1
            end
            
            for _ in 1:n_bank
                E_new = sample_fission_energy(xs, rng)
                push!(fission_bank, SourceSite(p.r, E_new, 1.0))
            end
            
            # Kill particle (fission is an absorption)
            p.alive = false
            
        elseif reaction == CAPTURE
            # Absorbed
            p.alive = false
            
        elseif reaction == ELASTIC
            # Elastic scattering
            nuc = get_dominant_fissile(mat)  # Use heaviest nuclide
            A = Float64(nuc.A)
            
            # Sample scattering angle in CM frame
            μ_cm = 2*rand(rng) - 1
            
            # Energy loss
            α = ((A - 1)/(A + 1))^2
            E_min = α * p.E
            p.E = E_min + rand(rng) * (p.E - E_min)  # Uniform in [αE, E]
            
            # Lab frame scattering cosine
            μ_lab = sample_elastic_cosine(A, p.E, rng)
            
            # Rotate direction
            p.u = rotate_direction(p.u, μ_lab, rng)
            
        elseif reaction == INELASTIC
            # Inelastic scattering - simplified treatment
            # Assume some energy loss and isotropic emission
            p.E *= 0.5 + 0.5*rand(rng)  # Lose 0-50% energy
            p.u = sample_isotropic_direction(rng)
        end
        
    else
        # Particle escapes
        
        # Track-length contribution up to boundary
        νΣf = get_νΣ_fission(mat, p.E)
        track_length_keff = p.weight * d_boundary * νΣf
        
        # Move to boundary and kill
        p.r = p.r + d_boundary * p.u
        p.alive = false
    end
    
    return (track_length_keff, collision_keff)
end

"""
    transport_particle!(p::Particle, geom::Geometry, mat::Material,
                       fission_bank::Vector{SourceSite}, rng::AbstractRNG,
                       k_eff::Float64) -> Tuple{Float64, Float64}

Transport a particle until it dies (absorption or escape).
Returns accumulated (track_length_tally, collision_tally).
"""
function transport_particle!(
    p::Particle,
    geom::Geometry,
    mat::Material,
    fission_bank::Vector{SourceSite},
    rng::AbstractRNG,
    k_eff::Float64
)
    tl_total = 0.0
    col_total = 0.0
    
    max_steps = 10000  # Prevent infinite loops
    step = 0
    
    while p.alive && step < max_steps
        tl, col = transport_step!(p, geom, mat, fission_bank, rng, k_eff)
        tl_total += tl
        col_total += col
        step += 1
    end
    
    return (tl_total, col_total)
end
