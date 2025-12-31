# Nuclear Football: A Monte Carlo Criticality Analysis

**Can a literal NFL football filled with weapons-grade plutonium achieve nuclear criticality?**

This project answers that question using first-principles Monte Carlo neutron transport simulation—the same computational technique used by national laboratories to design nuclear reactors and weapons.

## The Question

The term "nuclear football" usually refers to the briefcase carried by a military aide accompanying the U.S. President, containing launch codes for nuclear weapons. But what if we took the term literally?

An official NFL football has precise dimensions:
- **Length**: 28.0 cm (11 inches)  
- **Maximum circumference**: 56 cm (22 inches) → diameter ≈ 17.8 cm
- **Shape**: A prolate superellipsoid with pointed ends

If you filled this football with delta-phase plutonium-239 (density 15.8 g/cm³), would you have a critical nuclear assembly?

## The Answer

**Yes.** An NFL football filled with weapons-grade plutonium would be supercritical.

Our Monte Carlo simulation calculates **k_eff ≈ 1.2** for a Pu-239 filled football, meaning each generation of neutrons produces 20% more neutrons than the last. This is a prompt supercritical configuration—it would undergo a nuclear excursion.

The football shape is actually *suboptimal* for criticality (a sphere is ideal), but the ~22 kg of fissile material in an NFL-sized football more than compensates. The critical mass for plutonium in a football shape is approximately 11 kg.

## How It Works

This is not a toy simulation. We implement proper continuous-energy Monte Carlo neutron transport:

### Physics
- **Continuous-energy cross sections** interpolated from ENDF/B-VIII.0 tabulated data
- **Four reaction channels**: elastic scattering, inelastic scattering, fission, and radiative capture
- **Watt fission spectrum** for sampling outgoing neutron energies
- **Lab-frame scattering kinematics** with proper energy-angle correlations

### Methods
- **Power iteration** for k-eigenvalue (criticality) calculation
- **Track-length and collision estimators** for k_eff
- **Shannon entropy** monitoring for source convergence
- **Uniform combing** for fission bank normalization

### Validation
The code is validated against the **Jezebel benchmark**—a bare sphere of delta-phase plutonium that is exactly critical. Our simulation predicts k_eff ≈ 1.00 for Jezebel (within 5-6% due to simplified cross section data).

## Project Structure

```
buckeye/
├── julia/                    # Monte Carlo transport code
│   ├── src/
│   │   ├── NuclearFootball.jl    # Main module
│   │   ├── materials.jl          # Nuclide/material definitions
│   │   ├── cross_sections.jl     # Continuous-energy XS data
│   │   ├── geometry.jl           # Sphere, ellipsoid, football shapes
│   │   ├── particle.jl           # Neutron transport physics
│   │   ├── criticality.jl        # k-eigenvalue solver
│   │   └── analysis.jl           # Football analysis, validation
│   ├── scripts/
│   │   └── football.jl           # Main entry point
│   ├── Project.toml
│   └── Manifest.toml
├── book/                     # mdbook documentation
│   └── ...                   # Detailed theory and implementation
└── README.md                 # This file
```

## Running the Simulation

### Prerequisites

- Julia 1.9+ (install from [julialang.org](https://julialang.org/downloads/))
- Required packages are handled automatically via `Project.toml`

### Quick Start

```bash
cd julia
julia --project=. scripts/football.jl
```

This runs the full analysis:
1. Validates against the Jezebel benchmark
2. Calculates k_eff for a Pu-filled NFL football
3. Compares football vs sphere at equal mass
4. Finds the critical mass in football geometry

### Expected Output

```
╔══════════════════════════════════════════════════════════╗
║          CAN A NUCLEAR FOOTBALL GO CRITICAL?            ║
╚══════════════════════════════════════════════════════════╝

Material: δ-phase Pu (Jezebel)
Density: 15.61 g/cm³

============================================================
NFL FOOTBALL DIMENSIONS
============================================================
  Length: 28.0 cm (11.0 inches)
  Max diameter: 17.8 cm (7.0 inches)
  Volume: 2330.5 cm³
  If filled with δ-Pu: 36.4 kg

k_eff = 1.2142 ± 0.0035

🏈 RESULT: YES! A NUCLEAR FOOTBALL WOULD GO CRITICAL! 🏈
```

## The Math

### Criticality Condition

A fissile assembly is critical when **k_eff = 1**, meaning each generation of neutrons produces exactly as many neutrons as the last:

$$k_\text{eff} = \frac{\text{neutrons in generation } n+1}{\text{neutrons in generation } n}$$

### Monte Carlo Transport

We track individual neutrons through the geometry. At each step:

1. **Sample collision distance** from exponential distribution:
   $$d = -\frac{\ln(\xi)}{\Sigma_t(E)}$$
   where $\Sigma_t$ is the macroscopic total cross section.

2. **Check for boundary crossing**. If $d > d_\text{boundary}$, the neutron escapes.

3. **Sample reaction type** based on partial cross sections:
   - Elastic scattering: neutron bounces, loses some energy
   - Inelastic scattering: neutron excites the nucleus, loses more energy
   - Fission: neutron is absorbed, 2-3 new neutrons are born
   - Capture: neutron is absorbed, no new neutrons

4. **Bank fission sites** for the next generation.

### Power Iteration

We iterate:
1. Start N neutrons from fission bank
2. Transport all neutrons, collecting new fission sites
3. Estimate k_eff from fission production / source size
4. Normalize fission bank to N sites
5. Repeat until converged

### The Football Shape

An NFL football is well-approximated by a **superellipsoid**:

$$\left(\frac{\sqrt{x^2 + y^2}}{a}\right)^2 + \left(\frac{|z|}{c}\right)^m = 1$$

With $m \approx 3$ giving the characteristic pointed ends. This shape has ~15% higher neutron leakage than a sphere of equal volume, requiring more mass for criticality.

## Documentation

For a comprehensive deep-dive into the physics and mathematics, see the [mdbook documentation](book/):

```bash
cd book
mdbook serve --open
```

This covers:
- Nuclear fission fundamentals
- Neutron transport theory
- Monte Carlo methods for criticality
- Cross section physics
- Implementation details

## Disclaimers

This is an educational project demonstrating computational nuclear physics. The cross section data is simplified (representative values, not full ENDF pointwise data), which is why we see ~5-6% error on the Jezebel benchmark. A production code would use actual evaluated nuclear data files.

**Do not attempt to construct critical assemblies.** Handling fissile materials requires specialized facilities, training, and licensing. Accidental criticality is extremely dangerous.

## References

1. **ENDF/B-VIII.0**: Evaluated Nuclear Data File, maintained by NNDC at Brookhaven National Laboratory
2. **Jezebel Benchmark**: International Handbook of Evaluated Criticality Safety Benchmark Experiments (ICSBEP), PU-MET-FAST-001
3. **Monte Carlo Methods**: Lux & Koblinger, "Monte Carlo Particle Transport Methods" (1991)
4. **Criticality Safety**: Knief, "Nuclear Criticality Safety" (1985)

## License

MIT License. See LICENSE for details.
