"""
    Geometry definitions

Defines geometric shapes for Monte Carlo transport:
- Sphere (for validation)
- Ellipsoid (prolate spheroid)
- Football (realistic NFL football shape)

All geometries support:
- Point containment test
- Distance to boundary
- Surface normal
- Volume calculation
- Uniform point sampling
"""

#=============================================================================
# Abstract Geometry Interface
=============================================================================#

"""
    Geometry

Abstract type for all geometric shapes.
"""
abstract type Geometry end

"""
    point_inside(geom::Geometry, r::Vec3) -> Bool

Test if point r is inside the geometry.
"""
function point_inside end

"""
    distance_to_boundary(geom::Geometry, r::Vec3, u::Vec3) -> Float64

Distance from point r along direction u to the boundary.
Returns Inf if no intersection.
"""
function distance_to_boundary end

"""
    surface_normal(geom::Geometry, r::Vec3) -> Vec3

Outward unit normal at surface point r.
"""
function surface_normal end

"""
    volume(geom::Geometry) -> Float64

Volume of the geometry in cm³.
"""
function volume end

"""
    bounding_box(geom::Geometry) -> Tuple{Vec3, Vec3}

Axis-aligned bounding box: (min_corner, max_corner).
"""
function bounding_box end

"""
    sample_point(geom::Geometry, rng::AbstractRNG) -> Vec3

Sample a uniformly distributed random point inside the geometry.
"""
function sample_point end

#=============================================================================
# Sphere
=============================================================================#

"""
    Sphere

A sphere centered at the origin.
"""
struct Sphere <: Geometry
    radius::Float64  # cm
end

function point_inside(s::Sphere, r::Vec3)
    return norm(r) ≤ s.radius
end

function distance_to_boundary(s::Sphere, r::Vec3, u::Vec3)
    # Ray: P(t) = r + t*u
    # Sphere: |P|² = R²
    # |r + tu|² = R²
    # t² + 2(r·u)t + (|r|² - R²) = 0
    
    a = 1.0  # |u|² = 1 for unit vector
    b = 2.0 * dot(r, u)
    c = dot(r, r) - s.radius^2
    
    discriminant = b^2 - 4*a*c
    
    if discriminant < 0
        return Inf
    end
    
    sqrt_disc = sqrt(discriminant)
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    
    # We want the smallest positive t
    if t1 > 1e-10
        return t1
    elseif t2 > 1e-10
        return t2
    else
        return Inf
    end
end

function surface_normal(s::Sphere, r::Vec3)
    return normalize(r)
end

function volume(s::Sphere)
    return (4/3) * π * s.radius^3
end

function bounding_box(s::Sphere)
    R = s.radius
    return (Vec3(-R, -R, -R), Vec3(R, R, R))
end

function sample_point(s::Sphere, rng::AbstractRNG)
    # Rejection sampling in bounding box
    R = s.radius
    while true
        x = (2*rand(rng) - 1) * R
        y = (2*rand(rng) - 1) * R
        z = (2*rand(rng) - 1) * R
        r = Vec3(x, y, z)
        if norm(r) ≤ R
            return r
        end
    end
end

#=============================================================================
# Ellipsoid (Prolate Spheroid)
=============================================================================#

"""
    Ellipsoid

An ellipsoid centered at origin with semi-axes (a, b, c).
Surface equation: (x/a)² + (y/b)² + (z/c)² = 1

For a prolate spheroid (football-like): a = b < c (elongated along z)
"""
struct Ellipsoid <: Geometry
    a::Float64  # semi-axis in x direction (cm)
    b::Float64  # semi-axis in y direction (cm)  
    c::Float64  # semi-axis in z direction (cm)
end

"""
    Ellipsoid(equatorial::Float64, polar::Float64)

Create a prolate spheroid with given equatorial and polar radii.
Elongated along z-axis.
"""
function Ellipsoid(equatorial::Float64, polar::Float64)
    return Ellipsoid(equatorial, equatorial, polar)
end

function point_inside(e::Ellipsoid, r::Vec3)
    return (r[1]/e.a)^2 + (r[2]/e.b)^2 + (r[3]/e.c)^2 ≤ 1.0
end

function distance_to_boundary(e::Ellipsoid, r::Vec3, u::Vec3)
    # Transform to unit sphere coordinates
    # x' = x/a, y' = y/b, z' = z/c
    # Then solve ray-sphere intersection
    
    r_scaled = Vec3(r[1]/e.a, r[2]/e.b, r[3]/e.c)
    u_scaled = Vec3(u[1]/e.a, u[2]/e.b, u[3]/e.c)
    
    # |r' + t*u'|² = 1
    a_coef = dot(u_scaled, u_scaled)
    b_coef = 2.0 * dot(r_scaled, u_scaled)
    c_coef = dot(r_scaled, r_scaled) - 1.0
    
    discriminant = b_coef^2 - 4*a_coef*c_coef
    
    if discriminant < 0
        return Inf
    end
    
    sqrt_disc = sqrt(discriminant)
    t1 = (-b_coef - sqrt_disc) / (2*a_coef)
    t2 = (-b_coef + sqrt_disc) / (2*a_coef)
    
    if t1 > 1e-10
        return t1
    elseif t2 > 1e-10
        return t2
    else
        return Inf
    end
end

function surface_normal(e::Ellipsoid, r::Vec3)
    # Gradient of F(x,y,z) = (x/a)² + (y/b)² + (z/c)² - 1
    # ∇F = (2x/a², 2y/b², 2z/c²)
    grad = Vec3(2*r[1]/e.a^2, 2*r[2]/e.b^2, 2*r[3]/e.c^2)
    return normalize(grad)
end

function volume(e::Ellipsoid)
    return (4/3) * π * e.a * e.b * e.c
end

function bounding_box(e::Ellipsoid)
    return (Vec3(-e.a, -e.b, -e.c), Vec3(e.a, e.b, e.c))
end

function sample_point(e::Ellipsoid, rng::AbstractRNG)
    # Rejection sampling
    while true
        x = (2*rand(rng) - 1) * e.a
        y = (2*rand(rng) - 1) * e.b
        z = (2*rand(rng) - 1) * e.c
        r = Vec3(x, y, z)
        if point_inside(e, r)
            return r
        end
    end
end

#=============================================================================
# Football - Realistic NFL Football Shape
=============================================================================#

"""
    Football

A realistic American football shape using a superellipsoid.
The NFL specifies:
- Long axis: 11.0 - 11.25 inches (27.9 - 28.6 cm)
- Short circumference: 20.75 - 21.25 inches → diameter ~6.6-6.8 cm at center
- Pointed ends (not a true ellipsoid)

We model this as a superellipsoid:
    (|x|/a)^n + (|y|/a)^n + (|z|/c)^m = 1
    
With n=2 (circular cross-section) and m≈3-4 (pointed ends).
"""
struct Football <: Geometry
    half_length::Float64   # Half-length along z-axis (cm)
    max_radius::Float64    # Maximum radius at center (cm)
    pointiness::Float64    # Exponent m for z-axis (higher = pointier)
end

"""
    Football()

Create an official NFL football geometry.
- Length: 28.0 cm (11 inches)
- Max diameter: 17.8 cm (7 inches) - circumference ~56 cm
"""
function Football()
    return Football(14.0, 8.9, 3.0)  # half_length=14cm, radius=8.9cm, m=3
end

"""
    Football(length_cm::Float64, diameter_cm::Float64)

Create a football with specified overall length and maximum diameter.
"""
function Football(length_cm::Float64, diameter_cm::Float64)
    return Football(length_cm/2, diameter_cm/2, 3.0)
end

function _football_implicit(f::Football, r::Vec3)
    # Superellipsoid: (√(x² + y²)/a)² + (|z|/c)^m = 1
    # Returns <0 inside, =0 on surface, >0 outside
    rho = sqrt(r[1]^2 + r[2]^2)  # cylindrical radius
    z_term = (abs(r[3]) / f.half_length)^f.pointiness
    r_term = (rho / f.max_radius)^2
    return r_term + z_term - 1.0
end

function point_inside(f::Football, r::Vec3)
    return _football_implicit(f, r) ≤ 0.0
end

function distance_to_boundary(f::Football, r::Vec3, u::Vec3)
    # For superellipsoid, use Newton's method along ray
    # P(t) = r + t*u, find t where F(P(t)) = 0
    
    # First, quick bounding box check
    if abs(r[3]) > f.half_length * 2 && sign(r[3]) == sign(u[3])
        return Inf  # Heading away from football
    end
    
    # Use bisection/Newton hybrid for robustness
    # Start with bracketing
    t_min = 0.0
    t_max = 2.0 * (f.half_length + f.max_radius)  # Upper bound
    
    # Find bracket
    F0 = _football_implicit(f, r)
    
    if F0 > 0
        # Starting outside - find entry point
        t = 1e-6
        while t < t_max
            p = r + t * u
            if _football_implicit(f, p) < 0
                # Found inside point, bracket is [0, t]
                t_max = t
                break
            end
            t *= 2.0
        end
        if t >= t_max
            return Inf  # Never enters
        end
    end
    
    # Newton refinement
    t = (t_min + t_max) / 2
    for _ in 1:50
        p = r + t * u
        F = _football_implicit(f, p)
        
        if abs(F) < 1e-10
            return t > 1e-10 ? t : Inf
        end
        
        # Numerical gradient along ray
        dt = 1e-8
        F_plus = _football_implicit(f, r + (t + dt) * u)
        dF_dt = (F_plus - F) / dt
        
        if abs(dF_dt) < 1e-14
            break
        end
        
        t_new = t - F / dF_dt
        
        # Keep in bounds
        t_new = clamp(t_new, t_min, t_max)
        
        # Bisection fallback if Newton overshoots
        if t_new ≤ t_min || t_new ≥ t_max
            t_new = (t_min + t_max) / 2
        end
        
        t = t_new
    end
    
    return t > 1e-10 ? t : Inf
end

function surface_normal(f::Football, r::Vec3)
    # Gradient of F = (ρ/a)² + (|z|/c)^m - 1
    # ∂F/∂x = 2x/a²
    # ∂F/∂y = 2y/a²
    # ∂F/∂z = m * sign(z) * |z|^(m-1) / c^m
    
    rho = sqrt(r[1]^2 + r[2]^2)
    
    dF_dx = 2 * r[1] / f.max_radius^2
    dF_dy = 2 * r[2] / f.max_radius^2
    
    if abs(r[3]) > 1e-10
        dF_dz = f.pointiness * sign(r[3]) * abs(r[3])^(f.pointiness - 1) / f.half_length^f.pointiness
    else
        dF_dz = 0.0
    end
    
    grad = Vec3(dF_dx, dF_dy, dF_dz)
    n = norm(grad)
    return n > 1e-10 ? grad / n : Vec3(0, 0, 1)
end

function volume(f::Football)
    # Volume of superellipsoid: V = (2abc/mn) * B(1/m, 1/n + 1) * B(1/n, 1/n)
    # For our case with n=2, this simplifies
    # Approximate by numerical integration or use formula
    
    # Use Monte Carlo for accuracy
    N = 100000
    count = 0
    bb_min, bb_max = bounding_box(f)
    bb_vol = prod(bb_max - bb_min)
    
    rng = Random.default_rng()
    for _ in 1:N
        x = bb_min[1] + rand(rng) * (bb_max[1] - bb_min[1])
        y = bb_min[2] + rand(rng) * (bb_max[2] - bb_min[2])
        z = bb_min[3] + rand(rng) * (bb_max[3] - bb_min[3])
        if point_inside(f, Vec3(x, y, z))
            count += 1
        end
    end
    
    return bb_vol * count / N
end

function bounding_box(f::Football)
    return (Vec3(-f.max_radius, -f.max_radius, -f.half_length),
            Vec3(f.max_radius, f.max_radius, f.half_length))
end

function sample_point(f::Football, rng::AbstractRNG)
    # Rejection sampling in bounding box
    bb_min, bb_max = bounding_box(f)
    while true
        x = bb_min[1] + rand(rng) * (bb_max[1] - bb_min[1])
        y = bb_min[2] + rand(rng) * (bb_max[2] - bb_min[2])
        z = bb_min[3] + rand(rng) * (bb_max[3] - bb_min[3])
        r = Vec3(x, y, z)
        if point_inside(f, r)
            return r
        end
    end
end

#=============================================================================
# Utility Functions
=============================================================================#

"""
    sphere_from_mass(mat::Material, mass_kg::Float64) -> Sphere

Create a sphere with given material and mass.
"""
function sphere_from_mass(mat::Material, mass_kg::Float64)
    mass_g = mass_kg * 1000
    vol_cm3 = mass_g / mat.density
    radius = (3 * vol_cm3 / (4π))^(1/3)
    return Sphere(radius)
end

"""
    scale_geometry(geom::Geometry, factor::Float64) -> Geometry

Scale a geometry by a linear factor.
"""
function scale_geometry(s::Sphere, factor::Float64)
    return Sphere(s.radius * factor)
end

function scale_geometry(e::Ellipsoid, factor::Float64)
    return Ellipsoid(e.a * factor, e.b * factor, e.c * factor)
end

function scale_geometry(f::Football, factor::Float64)
    return Football(f.half_length * factor, f.max_radius * factor, f.pointiness)
end

"""
    geometry_with_volume(template::Geometry, target_volume::Float64) -> Geometry

Scale a geometry template to achieve target volume.
"""
function geometry_with_volume(template::Geometry, target_volume::Float64)
    current_volume = volume(template)
    scale_factor = (target_volume / current_volume)^(1/3)
    return scale_geometry(template, scale_factor)
end

"""
    geometry_with_mass(template::Geometry, mat::Material, mass_kg::Float64) -> Geometry

Scale a geometry to hold a specific mass of material.
"""
function geometry_with_mass(template::Geometry, mat::Material, mass_kg::Float64)
    mass_g = mass_kg * 1000
    target_volume = mass_g / mat.density
    return geometry_with_volume(template, target_volume)
end
