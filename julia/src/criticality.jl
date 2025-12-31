"""
    Criticality calculation

Monte Carlo k-eigenvalue (criticality) solver using power iteration.
"""

#=============================================================================
# Results Structure
=============================================================================#

"""
    KEffEstimates

k-effective estimates from different estimators.
"""
struct KEffEstimates
    track_length::Float64
    collision::Float64
    combined::Float64
end

"""
    CriticalityResult

Complete results from a criticality calculation.
"""
struct CriticalityResult
    k_eff::Float64
    k_std::Float64
    k_history::Vector{Float64}
    entropy_history::Vector{Float64}
    n_particles::Int
    n_generations::Int
    n_inactive::Int
    geometry::Geometry
    material::Material
end

#=============================================================================
# Shannon Entropy for Source Convergence
=============================================================================#

"""
    shannon_entropy(sites::Vector{SourceSite}, geom::Geometry, n_bins::Int=10) -> Float64

Calculate Shannon entropy of source distribution.
Higher entropy indicates more uniform distribution.
"""
function shannon_entropy(sites::Vector{SourceSite}, geom::Geometry, n_bins::Int=10)
    isempty(sites) && return 0.0
    
    # Get bounding box
    bb_min, bb_max = bounding_box(geom)
    
    # Create 3D mesh
    counts = zeros(Int, n_bins, n_bins, n_bins)
    total = 0
    
    for site in sites
        # Map position to bin
        ix = clamp(floor(Int, (site.r[1] - bb_min[1]) / (bb_max[1] - bb_min[1]) * n_bins) + 1, 1, n_bins)
        iy = clamp(floor(Int, (site.r[2] - bb_min[2]) / (bb_max[2] - bb_min[2]) * n_bins) + 1, 1, n_bins)
        iz = clamp(floor(Int, (site.r[3] - bb_min[3]) / (bb_max[3] - bb_min[3]) * n_bins) + 1, 1, n_bins)
        
        counts[ix, iy, iz] += 1
        total += 1
    end
    
    # Calculate entropy: H = -Σ p_i log₂(p_i)
    H = 0.0
    for c in counts
        if c > 0
            p = c / total
            H -= p * log2(p)
        end
    end
    
    return H
end

#=============================================================================
# Uniform Combing (Source Bank Normalization)
=============================================================================#

"""
    uniform_combing(fission_bank::Vector{SourceSite}, n_target::Int, rng::AbstractRNG) 
        -> Vector{SourceSite}

Select exactly n_target particles from fission bank using uniform combing.
This ensures reproducibility and avoids bias.
"""
function uniform_combing(fission_bank::Vector{SourceSite}, n_target::Int, rng::AbstractRNG)
    isempty(fission_bank) && error("Fission bank is empty - system is subcritical!")
    
    n_bank = length(fission_bank)
    
    if n_bank == n_target
        return copy(fission_bank)
    end
    
    # Calculate total weight
    total_weight = sum(site.weight for site in fission_bank)
    
    # Teeth spacing
    spacing = total_weight / n_target
    
    # Random starting offset
    offset = rand(rng) * spacing
    
    # Select sites
    new_bank = Vector{SourceSite}(undef, n_target)
    cumulative = 0.0
    tooth = offset
    site_idx = 1
    
    for i in 1:n_target
        # Find site containing this tooth
        while site_idx ≤ n_bank && cumulative + fission_bank[site_idx].weight < tooth
            cumulative += fission_bank[site_idx].weight
            site_idx += 1
        end
        
        if site_idx > n_bank
            site_idx = n_bank
        end
        
        # Copy site with unit weight
        site = fission_bank[site_idx]
        new_bank[i] = SourceSite(site.r, site.E, 1.0)
        
        tooth += spacing
    end
    
    return new_bank
end

#=============================================================================
# Main Criticality Solver
=============================================================================#

"""
    run_criticality(geom::Geometry, mat::Material;
                   n_particles::Int=10000,
                   n_inactive::Int=50,
                   n_active::Int=100,
                   show_progress::Bool=true) -> CriticalityResult

Run Monte Carlo k-eigenvalue calculation.

# Arguments
- `geom`: Geometry to simulate
- `mat`: Material filling the geometry
- `n_particles`: Neutrons per generation
- `n_inactive`: Generations to discard (source convergence)
- `n_active`: Generations to tally
- `show_progress`: Show progress bar
"""
function run_criticality(
    geom::Geometry,
    mat::Material;
    n_particles::Int=10000,
    n_inactive::Int=50,
    n_active::Int=100,
    show_progress::Bool=true
)
    rng = Random.default_rng()
    n_total = n_inactive + n_active
    
    # Get dominant fissile nuclide for initial energy sampling
    nuc = get_dominant_fissile(mat)
    xs = get_xs(nuc)
    
    # Initialize source bank with uniform distribution
    source_bank = Vector{SourceSite}(undef, n_particles)
    for i in 1:n_particles
        r = sample_point(geom, rng)
        E = sample_fission_energy(xs, rng)
        source_bank[i] = SourceSite(r, E, 1.0)
    end
    
    # Storage for results
    k_history = Float64[]
    entropy_history = Float64[]
    
    k_eff = 1.0  # Initial guess
    
    # Progress meter
    if show_progress
        prog = Progress(n_total, desc="Criticality: ", showspeed=true)
    end
    
    for gen in 1:n_total
        # New fission bank for this generation
        fission_bank = SourceSite[]
        
        # Tallies
        tl_tally = 0.0  # Track-length
        col_tally = 0.0 # Collision
        
        # Transport all particles
        for i in 1:n_particles
            site = source_bank[i]
            
            # Create particle from source site
            p = Particle(
                site.r,
                sample_isotropic_direction(rng),
                site.E,
                site.weight,
                true
            )
            
            # Transport until death
            tl, col = transport_particle!(p, geom, mat, fission_bank, rng, k_eff)
            tl_tally += tl
            col_tally += col
        end
        
        # Calculate k_eff estimates
        k_tl = tl_tally / n_particles
        k_col = col_tally / n_particles
        k_combined = 0.5 * (k_tl + k_col)
        
        # Use combined estimate for next generation
        k_eff = k_combined
        
        # Store history (only active generations for statistics)
        push!(k_history, k_eff)
        
        # Calculate entropy
        H = shannon_entropy(fission_bank, geom)
        push!(entropy_history, H)
        
        # Normalize fission bank for next generation
        if !isempty(fission_bank)
            source_bank = uniform_combing(fission_bank, n_particles, rng)
        else
            @warn "Empty fission bank at generation $gen - system may be subcritical"
            # Re-sample from geometry
            for i in 1:n_particles
                r = sample_point(geom, rng)
                E = sample_fission_energy(xs, rng)
                source_bank[i] = SourceSite(r, E, 1.0)
            end
        end
        
        if show_progress
            next!(prog, showvalues=[
                (:generation, gen),
                (:k_eff, @sprintf("%.5f", k_eff)),
                (:entropy, @sprintf("%.2f", H)),
                (:fission_sites, length(fission_bank))
            ])
        end
    end
    
    # Statistics from active generations only
    k_active = k_history[(n_inactive+1):end]
    k_mean = mean(k_active)
    k_std = std(k_active) / sqrt(n_active)  # Standard error
    
    return CriticalityResult(
        k_mean,
        k_std,
        k_history,
        entropy_history,
        n_particles,
        n_total,
        n_inactive,
        geom,
        mat
    )
end

#=============================================================================
# Critical Mass Search
=============================================================================#

"""
    find_critical_mass(template::Geometry, mat::Material;
                      mass_range::Tuple{Float64,Float64}=(1.0, 100.0),
                      tolerance::Float64=0.01,
                      n_particles::Int=5000,
                      n_inactive::Int=30,
                      n_active::Int=50) -> NamedTuple

Find the critical mass for a given geometry shape and material.

Returns a NamedTuple with:
- `mass_kg`: Critical mass in kg
- `k_eff`: k_eff at critical mass (should be ~1.0)
- `geometry`: The critical geometry
- `result`: Full CriticalityResult
"""
function find_critical_mass(
    template::Geometry,
    mat::Material;
    mass_range::Tuple{Float64,Float64}=(1.0, 100.0),
    tolerance::Float64=0.01,
    n_particles::Int=5000,
    n_inactive::Int=30,
    n_active::Int=50
)
    mass_lo, mass_hi = mass_range
    
    println("Searching for critical mass in range [$mass_lo, $mass_hi] kg")
    println("Geometry: $(typeof(template))")
    println("Material: $(mat.name)")
    println()
    
    # Bisection search
    iteration = 0
    max_iterations = 20
    
    while (mass_hi - mass_lo) / mass_lo > tolerance && iteration < max_iterations
        iteration += 1
        mass_mid = sqrt(mass_lo * mass_hi)  # Geometric mean
        
        # Create geometry at this mass
        geom = geometry_with_mass(template, mat, mass_mid)
        
        # Run criticality calculation
        result = run_criticality(
            geom, mat,
            n_particles=n_particles,
            n_inactive=n_inactive,
            n_active=n_active,
            show_progress=false
        )
        
        k = result.k_eff
        
        @printf("  Iteration %2d: mass = %6.2f kg, k_eff = %.4f ± %.4f\n",
                iteration, mass_mid, k, result.k_std)
        
        if k < 1.0
            mass_lo = mass_mid
        else
            mass_hi = mass_mid
        end
    end
    
    # Final calculation at converged mass
    final_mass = sqrt(mass_lo * mass_hi)
    final_geom = geometry_with_mass(template, mat, final_mass)
    
    println("\nFinal calculation at $(round(final_mass, digits=2)) kg...")
    final_result = run_criticality(
        final_geom, mat,
        n_particles=n_particles*2,
        n_inactive=n_inactive,
        n_active=n_active*2,
        show_progress=true
    )
    
    return (
        mass_kg = final_mass,
        k_eff = final_result.k_eff,
        k_std = final_result.k_std,
        geometry = final_geom,
        result = final_result
    )
end

#=============================================================================
# Quick k_eff Calculation
=============================================================================#

"""
    quick_k_eff(geom::Geometry, mat::Material; n_particles::Int=5000) -> Float64

Quick k_eff estimate with minimal generations.
Useful for parameter sweeps.
"""
function quick_k_eff(geom::Geometry, mat::Material; n_particles::Int=5000)
    result = run_criticality(
        geom, mat,
        n_particles=n_particles,
        n_inactive=20,
        n_active=30,
        show_progress=false
    )
    return result.k_eff
end
