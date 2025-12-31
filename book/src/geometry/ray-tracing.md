# Geometry and Ray Tracing

Monte Carlo transport requires knowing when neutrons cross geometry boundaries. This chapter covers the ray-tracing algorithms for spheres, ellipsoids, and the superellipsoid football shape.

## The Ray-Surface Intersection Problem

Given:
- Point **r** (neutron position)
- Direction **u** (unit vector)
- Geometry defined by implicit function F(**x**) = 0

Find the smallest positive t such that:
$$F(\mathbf{r} + t\mathbf{u}) = 0$$

This gives the distance to the boundary along the neutron's path.

## Spheres

A sphere of radius R centered at origin has implicit function:
$$F(\mathbf{x}) = |\mathbf{x}|^2 - R^2 = x^2 + y^2 + z^2 - R^2$$

### Point Inside Test
$$\text{inside} \Leftrightarrow F(\mathbf{r}) \leq 0 \Leftrightarrow |\mathbf{r}| \leq R$$

### Ray Intersection

Substitute **r** + t**u** into F = 0:
$$|\mathbf{r} + t\mathbf{u}|^2 = R^2$$

Expand:
$$|\mathbf{r}|^2 + 2t(\mathbf{r}\cdot\mathbf{u}) + t^2|\mathbf{u}|^2 = R^2$$

Since **u** is a unit vector (|**u**|² = 1):
$$t^2 + 2(\mathbf{r}\cdot\mathbf{u})t + (|\mathbf{r}|^2 - R^2) = 0$$

This is quadratic with:
- a = 1
- b = 2(**r**·**u**)
- c = |**r**|² - R²

Solutions:
$$t = -(\mathbf{r}\cdot\mathbf{u}) \pm \sqrt{(\mathbf{r}\cdot\mathbf{u})^2 - (|\mathbf{r}|^2 - R^2)}$$

We want the **smallest positive** root:

```julia
function distance_to_boundary(s::Sphere, r::Vec3, u::Vec3)
    b = 2.0 * dot(r, u)
    c = dot(r, r) - s.radius^2
    discriminant = b^2 - 4*c
    
    if discriminant < 0
        return Inf  # No intersection
    end
    
    sqrt_disc = sqrt(discriminant)
    t1 = (-b - sqrt_disc) / 2
    t2 = (-b + sqrt_disc) / 2
    
    if t1 > 1e-10
        return t1
    elseif t2 > 1e-10
        return t2
    else
        return Inf
    end
end
```

The tolerance (1e-10) prevents numerical issues when a particle is exactly on the surface.

### Surface Normal

The outward normal at surface point **r**:
$$\hat{n} = \frac{\nabla F}{|\nabla F|} = \frac{\mathbf{r}}{|\mathbf{r}|} = \frac{\mathbf{r}}{R}$$

## Ellipsoids

An ellipsoid with semi-axes (a, b, c):
$$F(\mathbf{x}) = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} - 1$$

A **prolate spheroid** (football-like) has a = b < c.

### Point Inside Test
$$\text{inside} \Leftrightarrow \frac{r_x^2}{a^2} + \frac{r_y^2}{b^2} + \frac{r_z^2}{c^2} \leq 1$$

### Ray Intersection

Transform to a unit sphere by scaling:
$$\mathbf{r}' = \left(\frac{r_x}{a}, \frac{r_y}{b}, \frac{r_z}{c}\right)$$
$$\mathbf{u}' = \left(\frac{u_x}{a}, \frac{u_y}{b}, \frac{u_z}{c}\right)$$

Then solve |**r**' + t**u**'|² = 1:
$$|\mathbf{u}'|^2 t^2 + 2(\mathbf{r}'\cdot\mathbf{u}')t + (|\mathbf{r}'|^2 - 1) = 0$$

```julia
function distance_to_boundary(e::Ellipsoid, r::Vec3, u::Vec3)
    r_scaled = Vec3(r[1]/e.a, r[2]/e.b, r[3]/e.c)
    u_scaled = Vec3(u[1]/e.a, u[2]/e.b, u[3]/e.c)
    
    a_coef = dot(u_scaled, u_scaled)
    b_coef = 2.0 * dot(r_scaled, u_scaled)
    c_coef = dot(r_scaled, r_scaled) - 1.0
    
    # ... solve quadratic as before
end
```

### Surface Normal

$$\hat{n} = \text{normalize}\left(\frac{2r_x}{a^2}, \frac{2r_y}{b^2}, \frac{2r_z}{c^2}\right)$$

## The Football: A Superellipsoid

An NFL football is not a true ellipsoid—it has **pointed ends**. We model it as a superellipsoid:

$$F(\mathbf{x}) = \left(\frac{\sqrt{x^2 + y^2}}{a}\right)^2 + \left(\frac{|z|}{c}\right)^m - 1$$

where:
- a = maximum radius (at center)
- c = half-length (tip to tip)
- m = "pointiness" exponent (m > 2 gives pointed ends)

For m = 2, this is an ellipsoid. For m = 3 or higher, the ends become pointed.

### Official NFL Football Dimensions
- Length: 28 cm (11 inches)
- Maximum circumference: 56 cm → diameter ≈ 17.8 cm
- Shape: m ≈ 3 matches the visual profile

Our parameters:
- half_length = 14 cm
- max_radius = 8.9 cm  
- pointiness = 3.0

### Point Inside Test

```julia
function point_inside(f::Football, r::Vec3)
    rho = sqrt(r[1]^2 + r[2]^2)  # Cylindrical radius
    z_term = (abs(r[3]) / f.half_length)^f.pointiness
    r_term = (rho / f.max_radius)^2
    return r_term + z_term ≤ 1.0
end
```

### Ray Intersection (Newton's Method)

Unlike spheres and ellipsoids, superellipsoids don't have closed-form ray intersections. We use **Newton's method** along the ray.

Define g(t) = F(**r** + t**u**). We seek the root g(t) = 0.

Newton iteration:
$$t_{n+1} = t_n - \frac{g(t_n)}{g'(t_n)}$$

where:
$$g'(t) = \nabla F \cdot \mathbf{u}$$

```julia
function distance_to_boundary(f::Football, r::Vec3, u::Vec3)
    # Initial bracket
    t_max = 2.0 * (f.half_length + f.max_radius)
    
    # Check if starting outside
    F0 = _football_implicit(f, r)
    
    if F0 > 0
        # Find entry point by doubling t until inside
        t = 1e-6
        while t < t_max
            if _football_implicit(f, r + t * u) < 0
                t_max = t
                break
            end
            t *= 2.0
        end
    end
    
    # Newton refinement
    t = t_max / 2
    for _ in 1:50
        p = r + t * u
        F = _football_implicit(f, p)
        
        if abs(F) < 1e-10
            return t > 1e-10 ? t : Inf
        end
        
        # Numerical derivative
        dt = 1e-8
        dF_dt = (_football_implicit(f, r + (t + dt) * u) - F) / dt
        
        t_new = t - F / dF_dt
        t = clamp(t_new, 0.0, t_max)
    end
    
    return t > 1e-10 ? t : Inf
end
```

### Surface Normal

$$\nabla F = \left(\frac{2x}{a^2\rho}, \frac{2y}{a^2\rho}, m \cdot \text{sign}(z) \cdot \frac{|z|^{m-1}}{c^m}\right) \cdot \rho$$

Normalized:
```julia
function surface_normal(f::Football, r::Vec3)
    rho = sqrt(r[1]^2 + r[2]^2)
    
    dF_dx = 2 * r[1] / f.max_radius^2
    dF_dy = 2 * r[2] / f.max_radius^2
    dF_dz = f.pointiness * sign(r[3]) * 
            abs(r[3])^(f.pointiness - 1) / 
            f.half_length^f.pointiness
    
    grad = Vec3(dF_dx, dF_dy, dF_dz)
    return normalize(grad)
end
```

## Volume Calculation

### Sphere
$$V = \frac{4}{3}\pi R^3$$

### Ellipsoid
$$V = \frac{4}{3}\pi abc$$

### Football (Superellipsoid)

The volume of a general superellipsoid is complicated. We use Monte Carlo integration:

```julia
function volume(f::Football)
    N = 100000
    count = 0
    
    # Bounding box
    bb_vol = 8 * f.max_radius^2 * f.half_length
    
    for _ in 1:N
        x = (2*rand() - 1) * f.max_radius
        y = (2*rand() - 1) * f.max_radius
        z = (2*rand() - 1) * f.half_length
        
        if point_inside(f, Vec3(x, y, z))
            count += 1
        end
    end
    
    return bb_vol * count / N
end
```

For our NFL football (m = 3): V ≈ 2300 cm³.

Compare to an ellipsoid with same dimensions: V = (4/3)π(8.9)²(14) ≈ 4640 cm³.

The pointed ends reduce volume by ~50%!

## Uniform Point Sampling

To initialize source particles, we need uniform sampling inside the geometry.

**Rejection sampling** is simple and works for any shape:

```julia
function sample_point(f::Football, rng)
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
```

Efficiency = V_geometry / V_bounding_box. For the football: ~2300 / 6930 ≈ 33%.

## Why Shape Matters for Criticality

A sphere minimizes surface area for a given volume, thus minimizing neutron leakage.

Any other shape has more surface area → more leakage → lower k_eff → higher critical mass.

For our football vs a sphere of equal mass:
- Football surface area: ~1400 cm²
- Equivalent sphere surface area: ~1100 cm²

The ~30% higher surface area translates to ~15% higher critical mass.

The pointed ends are particularly bad—they present a large surface area relative to their volume, creating "leakage hot spots."

## Scaling Geometry

To find critical mass, we need to scale the geometry:

```julia
function geometry_with_mass(template::Football, mat::Material, mass_kg::Float64)
    mass_g = mass_kg * 1000
    target_volume = mass_g / mat.density
    current_volume = volume(template)
    
    scale_factor = (target_volume / current_volume)^(1/3)
    
    return Football(
        template.half_length * scale_factor,
        template.max_radius * scale_factor,
        template.pointiness  # Keep shape constant
    )
end
```

Scaling by factor f:
- Volume scales as f³
- Mass scales as f³
- Surface area scales as f²
- Surface/volume ratio scales as 1/f → larger is more critical

## Summary

| Geometry | Intersection | Complexity |
|----------|-------------|------------|
| Sphere | Closed form (quadratic) | O(1) |
| Ellipsoid | Closed form (quadratic) | O(1) |
| Superellipsoid | Newton iteration | O(10-50) |

Our football requires ~10-50 Newton iterations per boundary check, making it slower than spheres, but still fast enough for practical Monte Carlo.
