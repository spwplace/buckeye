"""
    Analysis and comparison tools

Functions for analyzing nuclear football criticality.
"""

#=============================================================================
# Pretty Printing
=============================================================================#

"""
    print_results(result::CriticalityResult)

Print formatted criticality results.
"""
function print_results(result::CriticalityResult)
    println()
    println("=" ^ 60)
    println("CRITICALITY RESULTS")
    println("=" ^ 60)
    
    @printf("  k_eff = %.5f ± %.5f\n", result.k_eff, result.k_std)
    println()
    
    if result.k_eff > 1.0
        @printf("  STATUS: SUPERCRITICAL (k > 1)\n")
        excess = (result.k_eff - 1.0) * 100
        @printf("  Excess reactivity: %.2f%%\n", excess)
    elseif result.k_eff < 1.0
        @printf("  STATUS: SUBCRITICAL (k < 1)\n")
        deficit = (1.0 - result.k_eff) * 100
        @printf("  Reactivity deficit: %.2f%%\n", deficit)
    else
        @printf("  STATUS: CRITICAL (k ≈ 1)\n")
    end
    
    println()
    println("Simulation parameters:")
    @printf("  Particles/generation: %d\n", result.n_particles)
    @printf("  Inactive generations: %d\n", result.n_inactive)
    @printf("  Active generations:   %d\n", result.n_generations - result.n_inactive)
    
    println()
    println("Geometry: $(typeof(result.geometry))")
    @printf("  Volume: %.2f cm³\n", volume(result.geometry))
    
    println()
    println("Material: $(result.material.name)")
    @printf("  Density: %.2f g/cm³\n", result.material.density)
    
    mass_g = volume(result.geometry) * result.material.density
    @printf("  Total mass: %.2f kg\n", mass_g / 1000)
    
    println("=" ^ 60)
end

#=============================================================================
# Geometry Comparison
=============================================================================#

"""
    compare_geometries(mat::Material, mass_kg::Float64;
                      n_particles::Int=5000) -> Dict

Compare k_eff for different geometry shapes at the same mass.
"""
function compare_geometries(
    mat::Material,
    mass_kg::Float64;
    n_particles::Int=5000
)
    println("=" ^ 60)
    println("GEOMETRY COMPARISON")
    println("=" ^ 60)
    @printf("Material: %s\n", mat.name)
    @printf("Mass: %.2f kg\n", mass_kg)
    println()
    
    results = Dict{String, Any}()
    
    # Sphere (optimal)
    println("Testing SPHERE...")
    sphere = sphere_from_mass(mat, mass_kg)
    r_sphere = sphere.radius
    @printf("  Radius: %.2f cm\n", r_sphere)
    
    k_sphere = quick_k_eff(sphere, mat, n_particles=n_particles)
    results["sphere"] = (geometry=sphere, k_eff=k_sphere)
    @printf("  k_eff = %.4f\n\n", k_sphere)
    
    # Ellipsoid (2:1 aspect ratio)
    println("Testing ELLIPSOID (2:1 prolate)...")
    # Volume = (4/3)πabc with a=b, c=2a
    # V = (4/3)π * a² * 2a = (8/3)πa³
    vol = mass_kg * 1000 / mat.density
    a_ellipsoid = (vol * 3 / (8π))^(1/3)
    c_ellipsoid = 2 * a_ellipsoid
    ellipsoid = Ellipsoid(a_ellipsoid, c_ellipsoid)
    @printf("  Semi-axes: a=b=%.2f cm, c=%.2f cm\n", a_ellipsoid, c_ellipsoid)
    
    k_ellipsoid = quick_k_eff(ellipsoid, mat, n_particles=n_particles)
    results["ellipsoid_2to1"] = (geometry=ellipsoid, k_eff=k_ellipsoid)
    @printf("  k_eff = %.4f\n\n", k_ellipsoid)
    
    # Football (NFL dimensions scaled to mass)
    println("Testing FOOTBALL (NFL shape)...")
    football_template = Football()
    football = geometry_with_mass(football_template, mat, mass_kg)
    @printf("  Half-length: %.2f cm, Max radius: %.2f cm\n", 
            football.half_length, football.max_radius)
    
    k_football = quick_k_eff(football, mat, n_particles=n_particles)
    results["football"] = (geometry=football, k_eff=k_football)
    @printf("  k_eff = %.4f\n\n", k_football)
    
    # Summary
    println("-" ^ 60)
    println("SUMMARY")
    println("-" ^ 60)
    @printf("%-20s %10s %10s\n", "Geometry", "k_eff", "vs Sphere")
    println("-" ^ 60)
    
    for (name, data) in sort(collect(results), by=x->-x[2].k_eff)
        diff = (data.k_eff / k_sphere - 1) * 100
        @printf("%-20s %10.4f %+9.1f%%\n", name, data.k_eff, diff)
    end
    
    println("=" ^ 60)
    
    return results
end

#=============================================================================
# The Main Event: Football Analysis
=============================================================================#

"""
    analyze_football(; mass_kg::Float64=20.0, 
                     detailed::Bool=true) -> NamedTuple

The main analysis: Can a nuclear football achieve criticality?

Returns detailed results including:
- k_eff for a full Pu football
- Comparison with sphere at same mass
- Critical mass in football shape
"""
function analyze_football(;
    mass_kg::Float64=20.0,
    detailed::Bool=true,
    n_particles::Int=10000
)
    mat = δ_phase_plutonium
    
    println()
    println("╔" * "═"^58 * "╗")
    println("║" * " "^10 * "CAN A NUCLEAR FOOTBALL GO CRITICAL?" * " "^12 * "║")
    println("╚" * "═"^58 * "╝")
    println()
    
    println("Material: $(mat.name)")
    @printf("Density: %.2f g/cm³\n", mat.density)
    println()
    
    # First, what does an NFL football look like?
    nfl = Football()
    nfl_vol = volume(nfl)
    nfl_mass = nfl_vol * mat.density / 1000  # kg
    
    println("=" ^ 60)
    println("NFL FOOTBALL DIMENSIONS")
    println("=" ^ 60)
    @printf("  Length: %.1f cm (%.1f inches)\n", nfl.half_length * 2, nfl.half_length * 2 / 2.54)
    @printf("  Max diameter: %.1f cm (%.1f inches)\n", nfl.max_radius * 2, nfl.max_radius * 2 / 2.54)
    @printf("  Volume: %.1f cm³\n", nfl_vol)
    @printf("  If filled with δ-Pu: %.1f kg\n", nfl_mass)
    println()
    
    # Calculate k_eff for Pu-filled NFL football
    println("=" ^ 60)
    println("CRITICALITY OF PU-FILLED NFL FOOTBALL")
    println("=" ^ 60)
    println()
    
    println("Running Monte Carlo simulation...")
    println("(This may take a minute)")
    println()
    
    result_nfl = run_criticality(
        nfl, mat,
        n_particles=n_particles,
        n_inactive=50,
        n_active=100,
        show_progress=true
    )
    
    println()
    @printf("k_eff = %.4f ± %.4f\n", result_nfl.k_eff, result_nfl.k_std)
    println()
    
    if result_nfl.k_eff > 1.0
        println("🏈 RESULT: YES! A NUCLEAR FOOTBALL WOULD GO CRITICAL! 🏈")
        println()
        @printf("The system is SUPERCRITICAL by %.1f%%\n", (result_nfl.k_eff - 1) * 100)
    else
        println("Result: A standard NFL football filled with Pu-239 is SUBCRITICAL")
        println()
        @printf("Missing %.1f%% reactivity to reach critical\n", (1 - result_nfl.k_eff) * 100)
    end
    println()
    
    # Compare with sphere
    println("=" ^ 60)
    println("COMPARISON WITH OPTIMAL SPHERE")
    println("=" ^ 60)
    println()
    
    sphere = sphere_from_mass(mat, nfl_mass)
    @printf("Sphere with same mass (%.1f kg):\n", nfl_mass)
    @printf("  Radius: %.2f cm\n", sphere.radius)
    
    result_sphere = run_criticality(
        sphere, mat,
        n_particles=n_particles,
        n_inactive=50,
        n_active=100,
        show_progress=true
    )
    
    println()
    @printf("Sphere k_eff = %.4f ± %.4f\n", result_sphere.k_eff, result_sphere.k_std)
    
    shape_penalty = (result_sphere.k_eff / result_nfl.k_eff - 1) * 100
    @printf("\nShape penalty (football vs sphere): %.1f%%\n", shape_penalty)
    println()
    println("The football shape is less efficient due to higher neutron leakage")
    println("from the pointed ends.")
    
    if detailed
        # Find critical mass in football shape
        println()
        println("=" ^ 60)
        println("CRITICAL MASS IN FOOTBALL SHAPE")
        println("=" ^ 60)
        println()
        
        crit = find_critical_mass(
            Football(), mat,
            mass_range=(5.0, 50.0),
            tolerance=0.02,
            n_particles=n_particles÷2,
            n_inactive=30,
            n_active=50
        )
        
        println()
        @printf("Critical mass (football shape): %.2f kg\n", crit.mass_kg)
        
        # Compare with sphere critical mass
        crit_sphere = find_critical_mass(
            Sphere(1.0), mat,
            mass_range=(5.0, 50.0),
            tolerance=0.02,
            n_particles=n_particles÷2,
            n_inactive=30,
            n_active=50
        )
        
        println()
        @printf("Critical mass (sphere):         %.2f kg\n", crit_sphere.mass_kg)
        @printf("Football requires %.1f%% more material\n", 
                (crit.mass_kg / crit_sphere.mass_kg - 1) * 100)
    end
    
    # Final verdict
    println()
    println("╔" * "═"^58 * "╗")
    println("║" * " "^20 * "FINAL VERDICT" * " "^24 * "║")
    println("╠" * "═"^58 * "╣")
    
    if result_nfl.k_eff > 1.0
        println("║                                                          ║")
        println("║   ✓ YES, a plutonium-filled NFL football would be a     ║")
        println("║     supercritical nuclear device.                        ║")
        println("║                                                          ║")
        println("║   The pointed shape is suboptimal but the ~$(round(Int, nfl_mass)) kg of     ║")
        println("║   fissile material more than compensates.                ║")
        println("║                                                          ║")
    else
        println("║                                                          ║")
        println("║   ✗ A standard NFL football is too small to go critical ║")
        println("║     with delta-phase plutonium.                          ║")
        println("║                                                          ║")
    end
    
    println("╚" * "═"^58 * "╝")
    println()
    
    return (
        nfl_result = result_nfl,
        sphere_result = result_sphere,
        nfl_mass_kg = nfl_mass,
        critical_mass_football = detailed ? crit.mass_kg : nothing,
        critical_mass_sphere = detailed ? crit_sphere.mass_kg : nothing
    )
end

#=============================================================================
# Validation Against Jezebel Benchmark  
=============================================================================#

"""
    validate_jezebel() -> CriticalityResult

Validate against the Jezebel benchmark (bare Pu sphere).
Jezebel: 6.385 cm radius delta-phase Pu, should give k_eff ≈ 1.0
"""
function validate_jezebel()
    println("=" ^ 60)
    println("JEZEBEL BENCHMARK VALIDATION")
    println("=" ^ 60)
    println()
    println("Jezebel: Bare Pu-239 sphere (delta-phase with Ga)")
    println("Reference: Critical radius = 6.385 cm")
    println("Expected: k_eff ≈ 1.000")
    println()
    
    # Use Jezebel configuration
    jezebel = Sphere(6.385)  # cm
    mat = δ_phase_plutonium
    
    mass = volume(jezebel) * mat.density / 1000
    @printf("Mass: %.2f kg\n", mass)
    println()
    
    println("Running simulation...")
    result = run_criticality(
        jezebel, mat,
        n_particles=20000,
        n_inactive=50,
        n_active=150,
        show_progress=true
    )
    
    println()
    @printf("Calculated k_eff = %.5f ± %.5f\n", result.k_eff, result.k_std)
    @printf("Expected k_eff   = 1.00000\n")
    @printf("Difference       = %.2f%%\n", (result.k_eff - 1.0) * 100)
    println()
    
    if abs(result.k_eff - 1.0) < 0.05
        println("✓ VALIDATION PASSED (within 5%)")
    else
        println("✗ VALIDATION FAILED (>5% error)")
        println("  Note: Some error expected due to simplified cross sections")
    end
    
    println("=" ^ 60)
    
    return result
end
