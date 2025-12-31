"""
    Materials and nuclide definitions

Defines nuclear materials with their isotopic compositions and densities.
"""

#=============================================================================
# Nuclide Definition
=============================================================================# 

"""
    Nuclide

A nuclear isotope with its basic properties.
"""
struct Nuclide
    name::String
    Z::Int              # Atomic number
    A::Int              # Mass number
    atomic_mass::Float64  # Atomic mass in AMU
end

# Common nuclides
const Pu239 = Nuclide("Pu-239", 94, 239, 239.0521634)
const Pu240 = Nuclide("Pu-240", 94, 240, 240.0538135)
const Pu241 = Nuclide("Pu-241", 94, 241, 241.0568515)
const U235  = Nuclide("U-235",  92, 235, 235.0439299)
const U238  = Nuclide("U-238",  92, 238, 238.0507882)
const Ga69  = Nuclide("Ga-69",  31, 69,  68.9255736)
const Ga71  = Nuclide("Ga-71",  31, 71,  70.9247013)

#=============================================================================
# Material Definition  
=============================================================================#

"""
    Material

A material composed of one or more nuclides at specified atom fractions.
"""
struct Material
    name::String
    density::Float64                    # g/cm³
    composition::Vector{Tuple{Nuclide, Float64}}  # (nuclide, atom_fraction)
end

"""
    number_density(mat::Material, nuc::Nuclide) -> Float64

Get the number density of a specific nuclide in atoms/barn-cm.
"""
function number_density(mat::Material, nuc::Nuclide)
    # Find this nuclide's atom fraction
    atom_frac = 0.0
    for (n, frac) in mat.composition
        if n.name == nuc.name
            atom_frac = frac
            break
        end
    end
    atom_frac == 0.0 && return 0.0
    
    # Calculate average atomic mass
    avg_mass = sum(n.atomic_mass * f for (n, f) in mat.composition)
    
    # N = ρ * Nₐ / M * atom_fraction
    # Convert to atoms/barn-cm: divide by 1e24
    return mat.density * AVOGADRO / avg_mass * atom_frac / 1e24
end

"""
    total_number_density(mat::Material) -> Float64

Total atom number density in atoms/barn-cm.
"""
function total_number_density(mat::Material)
    avg_mass = sum(n.atomic_mass * f for (n, f) in mat.composition)
    return mat.density * AVOGADRO / avg_mass / 1e24
end

#=============================================================================
# Pre-defined Materials
=============================================================================#

"""
Delta-phase plutonium (Jezebel benchmark configuration).
Stabilized with ~1 wt% gallium. Density: 15.61 g/cm³
"""
const δ_phase_plutonium = Material(
    "δ-phase Pu (Jezebel)",
    15.61,  # g/cm³
    [
        (Pu239, 0.9545),  # 95.45 atom%
        (Pu240, 0.0045),  # 0.45 atom% (typical weapons grade)
        (Ga69,  0.0246),  # Natural Ga: 60.1% Ga-69
        (Ga71,  0.0164),  # Natural Ga: 39.9% Ga-71
    ]
)

"""
Alpha-phase plutonium (pure, unalloyed).
Higher density but brittle and pyrophoric. Density: 19.86 g/cm³
"""
const α_phase_plutonium = Material(
    "α-phase Pu (pure)",
    19.86,  # g/cm³
    [
        (Pu239, 1.0),
    ]
)

"""
Weapons-grade plutonium (typical composition).
< 7% Pu-240, delta-phase stabilized.
"""
const weapons_grade_pu = Material(
    "Weapons-grade Pu",
    15.8,  # g/cm³ (typical delta-phase)
    [
        (Pu239, 0.934),
        (Pu240, 0.060),
        (Pu241, 0.006),
    ]
)

"""
Highly-enriched uranium (93% U-235).
"""
const HEU = Material(
    "HEU (93%)",
    18.7,  # g/cm³
    [
        (U235, 0.93),
        (U238, 0.07),
    ]
)
