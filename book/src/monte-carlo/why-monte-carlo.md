# Why Monte Carlo?

Monte Carlo methods use random sampling to solve problems that would be intractable by other means. In neutron transport, we simulate individual neutrons—tracking their paths, collisions, and fates—to build up statistical estimates of quantities like k_eff.

## The Boltzmann Transport Equation

Neutron behavior is governed by the **Boltzmann transport equation**, a 7-dimensional integro-differential equation:

$$\frac{1}{v}\frac{\partial\psi}{\partial t} + \mathbf{\Omega} \cdot \nabla\psi + \Sigma_t\psi = \int_{4\pi}\int_0^\infty \Sigma_s(\mathbf{r}, E' \rightarrow E, \mathbf{\Omega}' \rightarrow \mathbf{\Omega})\psi' dE' d\Omega' + Q$$

where:
- $\psi(\mathbf{r}, E, \mathbf{\Omega}, t)$ is the angular neutron flux
- $\Sigma_t$ is the total macroscopic cross section  
- $\Sigma_s$ is the differential scattering cross section
- $Q$ is the external/fission source

The equation says: the rate of change of neutrons equals neutrons scattered in minus neutrons scattered out/absorbed plus sources.

Solving this equation analytically is possible only for highly idealized cases (infinite homogeneous media, simple geometries). Real problems require numerical methods.

## Deterministic vs Monte Carlo

### Deterministic Methods
Discretize the equation on a mesh of space, energy, and angle, then solve the resulting system of equations. Examples:
- **Diffusion theory**: Simplifies angular dependence
- **Discrete ordinates (Sn)**: Discretizes direction
- **Collision probability (CP)**: Tracks neutron flights between regions

Pros:
- Systematic convergence
- Low variance for integral quantities
- Fast for simple problems

Cons:
- Ray effects in Sn methods
- Geometric discretization errors
- Curse of dimensionality (7D problem → huge matrices)
- Struggles with complex geometry

### Monte Carlo
Simulate individual neutrons, tracking their random walks through the geometry. Tally quantities of interest and estimate statistics.

Pros:
- Exact geometry representation (no discretization)
- Continuous energy (no group averaging)
- Natural parallelization
- Handles any complexity

Cons:
- Statistical uncertainty (∝ 1/√N)
- Slow for deep penetration problems
- Variance can be large for local quantities
- Requires many histories for convergence

For our nuclear football—a complex 3D shape with continuous-energy physics—Monte Carlo is the natural choice.

## The Monte Carlo Philosophy

The fundamental insight is that we can **estimate integrals by random sampling**.

Consider estimating an integral:
$$I = \int_a^b f(x) dx$$

Instead of numerical quadrature, we:
1. Sample N random points $x_1, ..., x_N$ uniformly in [a, b]
2. Estimate: $\hat{I} = \frac{b-a}{N}\sum_{i=1}^N f(x_i)$

By the law of large numbers, $\hat{I} \rightarrow I$ as $N \rightarrow \infty$.

The standard error is:
$$\sigma_{\hat{I}} = \frac{b-a}{\sqrt{N}}\sigma_f$$

Key observation: the error decreases as $1/\sqrt{N}$, regardless of dimension. This is why Monte Carlo wins for high-dimensional problems.

## Analog vs Non-Analog Monte Carlo

### Analog Monte Carlo
Every simulated event corresponds directly to physical reality:
- Sample collision distance from exponential distribution
- Sample reaction type from physical probabilities
- Kill particle when absorbed or escaped

This is simple and unbiased, but can be inefficient (many neutrons wasted on uninteresting paths).

### Non-Analog (Variance Reduction)
Use weighted sampling to focus computation on important regions:
- **Implicit capture**: Don't kill on absorption; reduce weight instead
- **Weight windows**: Keep particle weights in a target range
- **Source biasing**: Preferentially sample important source regions
- **Geometry splitting**: Create copies of particles entering important regions

Our code uses mostly analog Monte Carlo with implicit capture, which is sufficient for criticality calculations.

## Random Number Generation

Monte Carlo requires a stream of random numbers. Modern codes use **pseudorandom number generators** (PRNGs) that are:
- Deterministic (given a seed)
- Statistically indistinguishable from true randomness
- Fast
- Long period (cycles before repeating)

Common choices:
- **Linear congruential generators**: x_{n+1} = (ax_n + c) mod m
- **Mersenne Twister**: Period 2^19937 - 1
- **Xorshift/Xoroshiro**: Fast, excellent statistical properties

Julia's default RNG (Xoshiro256++) is excellent for Monte Carlo work.

For reproducibility, we can set a seed:
```julia
using Random
Random.seed!(12345)
```

## History-Based Tracking

A Monte Carlo neutron "history" tracks one neutron from birth to death:

```
HISTORY START
├── Born at fission site (x₀, y₀, z₀), energy E₀
├── Travel distance d₁, collide
│   └── Elastic scatter, new direction and energy
├── Travel distance d₂, collide  
│   └── Inelastic scatter, lose energy
├── Travel distance d₃, collide
│   └── FISSION: create 3 fission sites for next generation
HISTORY END (absorbed by fission)
```

Each history contributes to tallies:
- Track-length tally: contribution = weight × distance × νΣf
- Collision tally: contribution = weight × νΣf/Σt at collision

After many histories, the central limit theorem gives us:
- Mean estimate of k
- Standard deviation of mean

## Parallel Monte Carlo

Monte Carlo is "embarrassingly parallel"—histories are independent, so we can:
- Run on multiple CPU cores
- Scale across cluster nodes
- Achieve near-linear speedup

The only synchronization needed is:
- Combining tallies at the end
- Managing the fission bank between generations (for criticality)

Our Julia implementation naturally benefits from Julia's threading model.

## Statistical Convergence

How many histories are enough? The standard error of the mean is:

$$\sigma_{\bar{x}} = \frac{\sigma}{\sqrt{N}}$$

To halve the uncertainty, we need 4× more histories. For k_eff:
- 10,000 histories/generation: σ ~ 0.01
- 100,000 histories/generation: σ ~ 0.003
- 1,000,000 histories/generation: σ ~ 0.001

For our football analysis, we use 10,000-20,000 histories per generation over 100-150 active generations, giving uncertainties around 0.003-0.005 in k_eff.

## What Makes Monte Carlo Work for Us

1. **Complex geometry**: The football's superellipsoid shape would be painful to mesh for deterministic methods
2. **Continuous energy**: We capture the full energy dependence of cross sections
3. **Modest precision needs**: We need k_eff to ~1%, which is achievable with reasonable runtime
4. **Clear physical interpretation**: Each simulated neutron corresponds to a real neutron path

Next, we'll see exactly how we track neutrons through the geometry.
