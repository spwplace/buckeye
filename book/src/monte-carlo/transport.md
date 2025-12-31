# Neutron Transport

This chapter details the algorithm for tracking a single neutron through a fissile assembly. This is the inner loop of our Monte Carlo code—executed millions of times per simulation.

## The Transport Loop

Each neutron is tracked from birth until death (absorption or escape). The algorithm:

```
function transport_neutron(particle, geometry, material):
    while particle.alive:
        # 1. Calculate distance to next collision
        Σ_total = get_macroscopic_cross_section(material, particle.energy)
        d_collision = sample_exponential(1/Σ_total)
        
        # 2. Calculate distance to boundary
        d_boundary = distance_to_surface(geometry, particle.position, particle.direction)
        
        # 3. Compare distances
        if d_collision < d_boundary:
            # Collision occurs inside geometry
            move_particle(particle, d_collision)
            process_collision(particle, material, fission_bank)
        else:
            # Particle escapes
            move_particle(particle, d_boundary)
            particle.alive = false  # Leakage
```

Let's examine each step in detail.

## Step 1: Sampling Collision Distance

A neutron traveling through matter has a constant probability per unit length of colliding. This leads to an exponential distribution of collision distances.

### The Physics

If Σ is the macroscopic total cross section (probability of interaction per cm), then the probability of traveling distance d without interaction is:

$$P(\text{no collision in } d) = e^{-\Sigma d}$$

The probability density for first collision at distance d is:

$$p(d) = \Sigma e^{-\Sigma d}$$

This is an exponential distribution with mean λ = 1/Σ (the mean free path).

### Sampling from Exponential

To sample from \\(p(d) = \Sigma e^{-\Sigma d}\\), we use the **inverse transform method**:

1. Sample ξ uniformly from [0, 1]
2. Set d = F⁻¹(ξ) where F is the CDF

The CDF is:
$$F(d) = \int_0^d \Sigma e^{-\Sigma d'} dd' = 1 - e^{-\Sigma d}$$

Inverting:
$$d = -\frac{\ln(1-\xi)}{\Sigma} = -\frac{\ln(\xi)}{\Sigma}$$

The last equality holds because ξ and 1-ξ have the same distribution.

### Implementation

```julia
function sample_collision_distance(Σ_total::Float64, rng::AbstractRNG)
    ξ = rand(rng)
    while ξ == 0.0  # Avoid log(0)
        ξ = rand(rng)
    end
    return -log(ξ) / Σ_total
end
```

For Pu-239 at 1 MeV with Σ_total ≈ 0.25 cm⁻¹, the mean free path is ~4 cm.

## Step 2: Distance to Boundary

We must determine if the neutron will hit the geometry boundary before its sampled collision point.

### Ray-Surface Intersection

Given particle position **r** and direction **u** (unit vector), find the smallest positive t where:

$$\mathbf{r} + t\mathbf{u} \in \partial\text{Geometry}$$

For a sphere of radius R:
$$|\mathbf{r} + t\mathbf{u}|^2 = R^2$$

Expanding:
$$t^2 + 2(\mathbf{r}\cdot\mathbf{u})t + (|\mathbf{r}|^2 - R^2) = 0$$

This is a quadratic with solutions:
$$t = -(\mathbf{r}\cdot\mathbf{u}) \pm \sqrt{(\mathbf{r}\cdot\mathbf{u})^2 - |\mathbf{r}|^2 + R^2}$$

We take the smallest positive root.

For the football (superellipsoid), we use Newton's method along the ray—see the Geometry chapter.

## Step 3: Process Collision

When a collision occurs, we must determine what happens: scatter, absorb, or fission.

### Sampling Reaction Type

The probability of each reaction is proportional to its macroscopic cross section:

$$P(\text{elastic}) = \frac{\Sigma_\text{el}}{\Sigma_t}, \quad P(\text{inelastic}) = \frac{\Sigma_\text{in}}{\Sigma_t}, \quad \text{etc.}$$

To sample:
1. Generate ξ ∈ [0, 1]
2. Compare to cumulative probabilities

```julia
function sample_reaction(material, E, rng)
    Σ_total = get_Σ_total(material, E)
    Σ_el = get_Σ_elastic(material, E)
    Σ_in = get_Σ_inelastic(material, E)
    Σ_f = get_Σ_fission(material, E)
    
    ξ = rand(rng) * Σ_total
    
    if ξ < Σ_el
        return ELASTIC
    elseif ξ < Σ_el + Σ_in
        return INELASTIC
    elseif ξ < Σ_el + Σ_in + Σ_f
        return FISSION
    else
        return CAPTURE
    end
end
```

## Elastic Scattering

In elastic scattering, the neutron bounces off the nucleus. We need to determine:
1. The new direction
2. The new energy

### Scattering Angle

In the center-of-mass (CM) frame, elastic scattering off heavy nuclei is approximately isotropic:
$$P(\mu_{cm}) = \frac{1}{2}, \quad \mu_{cm} \in [-1, 1]$$

Sample: \\(\mu_{cm} = 2\xi - 1\\)

### Energy Loss

The lab-frame energy after collision depends on the CM scattering angle:

$$E' = E \cdot \frac{A^2 + 2A\mu_{cm} + 1}{(A+1)^2}$$

For Pu-239 (A = 239), even backscattering (μ = -1) only reduces energy to 0.983E. Heavy nuclei are poor moderators.

### Direction Rotation

To rotate the direction vector by scattering angle θ (where cos θ = μ):

1. Sample azimuthal angle φ uniformly in [0, 2π]
2. Rotate **u** by angle θ around an axis perpendicular to **u**

The rotation formula:
$$\mathbf{u}' = \mu\mathbf{u} + \sqrt{1-\mu^2}\left[\frac{u_x u_z \cos\phi - u_y \sin\phi}{\sqrt{1-u_z^2}}, \frac{u_y u_z \cos\phi + u_x \sin\phi}{\sqrt{1-u_z^2}}, -\sqrt{1-u_z^2}\cos\phi\right]$$

Special case when **u** is nearly parallel to z-axis:
$$\mathbf{u}' = [\sqrt{1-\mu^2}\cos\phi, \sqrt{1-\mu^2}\sin\phi, \text{sign}(u_z)\mu]$$

## Inelastic Scattering

Inelastic scattering excites the nucleus, so the neutron loses more energy than kinematics alone would predict. We make a simplified treatment:
- Energy loss: E' = E × (0.5 + 0.5ξ) (lose 0-50%)
- Direction: Isotropic in lab frame

A production code would use proper level densities and angular distributions from nuclear data.

## Fission

The most important reaction for criticality. When fission occurs:

1. The incident neutron is absorbed (remove from simulation)
2. New neutrons are born at the fission site
3. Store these for the next generation

### Number of Neutrons

The average number is ν̄(E). We want exactly this many on average, so we use stochastic rounding:

```julia
expected = ν̄ * weight / k_eff
n_emit = floor(Int, expected)
if rand() < (expected - n_emit)
    n_emit += 1
end
```

The division by k_eff normalizes the fission source. If k > 1, fewer neutrons are banked per fission to prevent exponential growth of the bank.

### Fission Neutron Energy

Each fission neutron is sampled from the Watt spectrum:
$$\chi(E) \propto e^{-E/a} \sinh\sqrt{bE}$$

See the Cross Sections chapter for the sampling algorithm.

### Fission Neutron Direction

Fission neutrons are emitted isotropically:
$$\mathbf{u} = [\sqrt{1-\mu^2}\cos\phi, \sqrt{1-\mu^2}\sin\phi, \mu]$$

where μ = 2ξ₁ - 1 and φ = 2πξ₂.

## Capture

Radiative capture absorbs the neutron without producing new neutrons. Simply mark the particle as dead.

In non-analog Monte Carlo, we might use **implicit capture**—reduce the particle's weight by σ_c/σ_t at each collision instead of killing it. This reduces variance but our analog approach is simpler.

## Track-Length and Collision Tallies

As neutrons move, we accumulate contributions to the k_eff estimate.

### Track-Length Estimator

Every path segment contributes to k:
$$k_{TL} = \sum_{\text{segments}} w \cdot d \cdot \nu\Sigma_f(E)$$

where w is weight, d is distance traveled, and νΣf is the neutron production cross section.

Physical interpretation: the expected number of fission neutrons produced along this path.

### Collision Estimator

Every collision contributes:
$$k_{col} = \sum_{\text{collisions}} w \cdot \frac{\nu\Sigma_f(E)}{\Sigma_t(E)}$$

Physical interpretation: at each collision, the probability of fission times the number of neutrons produced.

### Combined Estimator

We average the two for reduced variance:
$$k_{eff} = \frac{1}{2}(k_{TL} + k_{col})$$

## Putting It Together

Here's the complete transport function from our code:

```julia
function transport_step!(p::Particle, geom::Geometry, mat::Material, 
                        fission_bank::Vector{SourceSite}, rng, k_eff)
    # Get cross section at current energy
    Σ_total = get_Σ_total(mat, p.E)
    
    # Sample collision distance
    d_collision = sample_collision_distance(Σ_total, rng)
    
    # Distance to boundary
    d_boundary = distance_to_boundary(geom, p.r, p.u)
    
    if d_collision < d_boundary
        # Collision inside geometry
        νΣf = get_νΣ_fission(mat, p.E)
        track_length_keff = p.weight * d_collision * νΣf
        
        # Move to collision site
        p.r = p.r + d_collision * p.u
        
        # Sample and process reaction
        reaction = sample_reaction(mat, p.E, rng)
        collision_keff = p.weight * νΣf / Σ_total
        
        if reaction == FISSION
            # Bank fission neutrons
            ...
            p.alive = false
        elseif reaction == CAPTURE
            p.alive = false
        elseif reaction == ELASTIC
            # Update energy and direction
            ...
        elseif reaction == INELASTIC
            # Update energy and direction
            ...
        end
    else
        # Escape
        νΣf = get_νΣ_fission(mat, p.E)
        track_length_keff = p.weight * d_boundary * νΣf
        
        p.r = p.r + d_boundary * p.u
        p.alive = false
    end
    
    return (track_length_keff, collision_keff)
end
```

A neutron typically undergoes 5-20 collisions before being absorbed or escaping, depending on system size and materials.

## Summary

| Step | Sampling Method | Physical Basis |
|------|----------------|----------------|
| Collision distance | Exponential | Constant collision probability/length |
| Reaction type | Discrete | Cross section ratios |
| Scatter angle | Uniform in CM | Isotropic in CM for heavy nuclei |
| Fission neutron count | Stochastic rounding | Preserves expected value |
| Fission energy | Watt spectrum | Empirical fission physics |
| Fission direction | Isotropic | No preferred direction |

Next, we'll see how to combine many neutron histories into a k_eff calculation.
