#!/usr/bin/env julia
"""
    Nuclear Football Analysis

The main event: Can a nuclear football achieve criticality?

Run with:
    julia --project=julia julia/scripts/football.jl
"""

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using NuclearFootball

function main()
    println()
    println("╔" * "═"^58 * "╗")
    println("║" * " "^8 * "NUCLEAR FOOTBALL CRITICALITY ANALYSIS" * " "^11 * "║")
    println("║" * " "^58 * "║")
    println("║" * " "^4 * "Monte Carlo Neutron Transport Simulation" * " "^13 * "║")
    println("╚" * "═"^58 * "╝")
    println()
    
    # First, validate against Jezebel benchmark
    println("Step 1: Validating code against Jezebel benchmark...")
    println()
    jezebel_result = validate_jezebel()
    println()
    
    # Then answer THE question
    println("Step 2: The main event...")
    println()
    
    results = analyze_football(
        detailed=true,
        n_particles=10000
    )
    
    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
