# The k-Eigenvalue Problem

Finding k_eff is an eigenvalue problem: we seek the multiplication factor k and the corresponding neutron flux distribution φ such that the system is in steady state. Monte Carlo solves this by simulating generation after generation of neutrons.

## The Mathematical Problem

The time-independent transport equation with fission is:

$$\hat{L}\phi = \frac{1}{k}\hat{F}\phi$$

where:
- \\(\hat{L}\\) is the net loss operator (streaming + absorption - scatter)
- \\(\hat{F}\\) is the fission production operator
- \\(k\\) is the eigenvalue (multiplication factor)
- \\(\phi\\) is the flux (eigenvector)

In integral form:
$$\phi(\mathbf{r}, E, \mathbf{\Omega}) = \int K(\mathbf{r}'\rightarrow\mathbf{r}, E'\rightarrow E, \mathbf{\Omega}'\rightarrow\mathbf{\Omega})\phi(\mathbf{r}', E', \mathbf{\Omega}')d\mathbf{r}'dE'd\mathbf{\Omega}' + \frac{1}{k}F$$

The fission source F depends on φ, making this a nonlinear problem.

## Power Iteration

The standard solution method is **power iteration** (also called the source iteration method):

1. **Initialize**: Guess a fission source distribution S₀(r)
2. **Transport**: Given source Sₙ, calculate flux φₙ
3. **Fission production**: Calculate new source Sₙ₊₁ = Fφₙ
4. **Estimate k**: \\(k_n = \frac{\int S_{n+1}}{\int S_n}\\)
5. **Normalize**: Sₙ₊₁ ← Sₙ₊₁/kₙ
6. **Converge**: Repeat until kₙ and Sₙ converge

Power iteration converges to the **fundamental mode**—the eigenvalue with largest magnitude (k_eff) and its corresponding eigenfunction.

### Convergence

The convergence rate depends on the **dominance ratio**:
$$\rho = \frac{k_1}{k_0}$$

where k₀ is the fundamental eigenvalue and k₁ is the next largest. For most reactors, ρ > 0.9, meaning slow convergence.

The error in kₙ after n iterations is:
$$k_n - k_{eff} \propto \rho^n$$

## Monte Carlo Implementation

In Monte Carlo, power iteration becomes generation-by-generation simulation:

```
POWER ITERATION (Monte Carlo):

Generation 0:
  - Create N particles from uniform fission source
  - Transport all particles
  - Collect fission sites → Fission Bank (M sites)
  - k₀ = M / N

Generation 1:
  - Sample N particles from Fission Bank
  - Transport all particles
  - Collect fission sites → New Bank (M₁ sites)  
  - k₁ = M₁ / N

...continue for many generations...
```

### Source Convergence

The initial source guess is typically uniform random within the geometry. This must converge to the true fission distribution before we can trust k_eff estimates.

We run **inactive generations** (typically 50-100) to allow source convergence, then **active generations** (100-200) to accumulate statistics.

### The Fission Bank

Each generation produces a variable number of fission sites. We must normalize to keep exactly N particles per generation.

**Uniform combing** (also called systematic resampling) selects exactly N particles from the bank while preserving the spatial distribution:

```julia
function uniform_combing(fission_bank, N_target, rng)
    total_weight = sum(site.weight for site in fission_bank)
    spacing = total_weight / N_target
    offset = rand(rng) * spacing
    
    new_bank = []
    cumulative = 0.0
    tooth = offset
    idx = 1
    
    for i in 1:N_target
        while idx ≤ length(fission_bank) && 
              cumulative + fission_bank[idx].weight < tooth
            cumulative += fission_bank[idx].weight
            idx += 1
        end
        push!(new_bank, SourceSite(fission_bank[idx].r, ...))
        tooth += spacing
    end
    
    return new_bank
end
```

This is more statistically efficient than simple random sampling.

## Shannon Entropy

How do we know when the source has converged? We monitor the **Shannon entropy** of the spatial distribution:

$$H = -\sum_i p_i \log_2 p_i$$

where pᵢ is the fraction of source particles in spatial bin i.

Properties:
- H = 0 when all particles are in one bin (maximally concentrated)
- H = log₂(N_bins) when uniformly distributed (maximally spread)

We mesh the geometry and count source particles in each cell. When H stabilizes, the source is converged.

```julia
function shannon_entropy(sites, geometry, n_bins=10)
    # Create 3D mesh
    counts = zeros(Int, n_bins, n_bins, n_bins)
    
    for site in sites
        # Map position to bin
        ix, iy, iz = position_to_bin(site.r, geometry, n_bins)
        counts[ix, iy, iz] += 1
    end
    
    # Calculate entropy
    total = sum(counts)
    H = 0.0
    for c in counts
        if c > 0
            p = c / total
            H -= p * log2(p)
        end
    end
    
    return H
end
```

## k_eff Estimators

We collect two estimates each generation:

### Track-Length Estimator
Sum over all path segments:
$$k_{TL} = \frac{1}{N}\sum_{\text{histories}}\sum_{\text{segments}} w \cdot d \cdot \nu\Sigma_f(E)$$

This accumulates continuously as particles move.

### Collision Estimator
Sum over all collisions:
$$k_{col} = \frac{1}{N}\sum_{\text{histories}}\sum_{\text{collisions}} w \cdot \frac{\nu\Sigma_f(E)}{\Sigma_t(E)}$$

### Combined Estimator
Average the two:
$$k_{eff} = \frac{1}{2}(k_{TL} + k_{col})$$

The combined estimator typically has lower variance than either alone.

## Statistical Analysis

After discarding inactive generations, we have a sequence of k estimates: {k₁, k₂, ..., kₙ}.

### Mean
$$\bar{k} = \frac{1}{n}\sum_{i=1}^n k_i$$

### Standard Error
$$\sigma_{\bar{k}} = \frac{\sigma_k}{\sqrt{n}}$$

where σ_k is the standard deviation of the kᵢ values.

For k_eff ≈ 1 with 10,000 particles over 100 generations:
- σ_k ≈ 0.03 (single generation)
- σ_k̄ ≈ 0.003 (mean over generations)

### Autocorrelation

Successive generations are correlated (they share the same fission source). The true standard error is:
$$\sigma_{\bar{k},\text{true}} = \sigma_{\bar{k}} \sqrt{1 + 2\sum_{\tau=1}^{\infty} \rho(\tau)}$$

where ρ(τ) is the autocorrelation at lag τ. For strongly correlated sources, this can increase uncertainty by 50-100%.

Production codes account for this; our simple code ignores it (acceptably for educational purposes).

## Our Implementation

The complete criticality solver:

```julia
function run_criticality(geom, mat; 
                        n_particles=10000,
                        n_inactive=50,
                        n_active=100)
    
    # Initialize source bank uniformly
    source_bank = [random_source_site(geom) for _ in 1:n_particles]
    
    k_history = Float64[]
    entropy_history = Float64[]
    k_eff = 1.0  # Initial guess
    
    for gen in 1:(n_inactive + n_active)
        fission_bank = SourceSite[]
        tl_tally = 0.0
        col_tally = 0.0
        
        # Transport all particles
        for site in source_bank
            particle = create_particle(site)
            while particle.alive
                tl, col = transport_step!(particle, geom, mat, 
                                         fission_bank, rng, k_eff)
                tl_tally += tl
                col_tally += col
            end
        end
        
        # Estimate k_eff
        k_tl = tl_tally / n_particles
        k_col = col_tally / n_particles
        k_eff = 0.5 * (k_tl + k_col)
        
        push!(k_history, k_eff)
        push!(entropy_history, shannon_entropy(fission_bank, geom))
        
        # Prepare next generation
        source_bank = uniform_combing(fission_bank, n_particles, rng)
    end
    
    # Statistics from active generations
    k_active = k_history[(n_inactive+1):end]
    k_mean = mean(k_active)
    k_std = std(k_active) / sqrt(n_active)
    
    return (k_eff=k_mean, k_std=k_std, history=k_history)
end
```

## Typical Output

```
Generation  k_eff   Entropy   Fission Sites
------------------------------------------
    1       0.987    4.21         9870
    2       1.023    4.35        10230
   ...
   50       1.182    5.89        11820  ← Source converged
   51       1.195    5.91        11950  ← Active generations start
  ...
  150       1.215    5.88        12150

Final: k_eff = 1.2034 ± 0.0041
```

The entropy stabilizes around generation 30-50, indicating source convergence. We then average over 100 active generations to get our final k_eff estimate.

## Summary

| Component | Purpose |
|-----------|---------|
| Power iteration | Converge to fundamental mode |
| Inactive generations | Allow source convergence |
| Active generations | Accumulate statistics |
| Shannon entropy | Monitor source convergence |
| Track-length estimator | Contribution from all paths |
| Collision estimator | Contribution at collisions |
| Uniform combing | Normalize fission bank |
| Standard error | Uncertainty quantification |

Our football simulation runs 50 inactive + 100 active generations with 10,000-20,000 particles, giving k_eff uncertainties of ~0.3-0.5%.
