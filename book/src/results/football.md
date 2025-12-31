# The Nuclear Football

The main event: analyzing criticality of an NFL football filled with weapons-grade plutonium.

## The Setup

### NFL Football Specifications

An official NFL football must meet these requirements (Rule 2, Section 1):
- **Long axis**: 11.0 to 11.25 inches (27.9 to 28.6 cm)
- **Long circumference**: 28.0 to 28.5 inches
- **Short circumference**: 21.0 to 21.25 inches → **diameter ≈ 17.8 cm**

We model this as a superellipsoid:
- Half-length: 14.0 cm
- Maximum radius: 8.9 cm
- Pointiness exponent: 3.0

### Material

We fill the football with **δ-phase plutonium-239**:
- Density: 15.61 g/cm³
- Composition: 95.45% Pu-239, 0.45% Pu-240, 4.1% Ga

This is the same material used in the Jezebel benchmark—stabilized with gallium to maintain the ductile delta phase.

### Calculated Properties

| Property | Value |
|----------|-------|
| Football volume | 2330 cm³ |
| Plutonium mass | 36.4 kg |
| Surface area | ~1400 cm² |

For comparison, a sphere of equal mass would have:
- Radius: 8.0 cm
- Volume: 2140 cm³
- Surface area: ~800 cm²

The football has 75% more surface area than an equal-mass sphere.

## The Simulation

```julia
# Create NFL football geometry
nfl = Football()  # 28 cm × 17.8 cm, m=3

# Run high-statistics calculation
result = run_criticality(
    nfl, δ_phase_plutonium,
    n_particles = 10000,
    n_inactive = 50,
    n_active = 100,
    show_progress = true
)
```

### Output

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

============================================================
CRITICALITY OF PU-FILLED NFL FOOTBALL
============================================================

Running Monte Carlo simulation...
Criticality: 100%|████████████████████| Time: 0:01:52

k_eff = 1.2142 ± 0.0038

🏈 RESULT: YES! A NUCLEAR FOOTBALL WOULD GO CRITICAL! 🏈

The system is SUPERCRITICAL by 21.4%
```

## Results

### k-effective

$$k_{eff} = 1.214 \pm 0.004$$

This is **deeply supercritical**. The football would undergo a prompt critical excursion.

### Reactivity

$$\rho = \frac{k-1}{k} = \frac{0.214}{1.214} = 17.6\%$$

In "dollars" (where \\(1 = β = 0.0022\\) for Pu):
$$\rho = \frac{0.176}{0.0022} = 80 \text{ dollars}$$

The system is **80 dollars supercritical**—far into the prompt supercritical regime.

### Time Scale

The neutron generation time in fast plutonium systems is ~10 ns. The e-folding time:
$$T = \frac{\ell}{k-1} = \frac{10^{-8}}{0.214} = 47 \text{ ns}$$

The neutron population would double every 47 nanoseconds. In 1 microsecond (~20 doublings), the population increases by a factor of ~10⁶.

## Comparison with Sphere

We also calculated k_eff for a sphere of equal mass:

```
============================================================
COMPARISON WITH OPTIMAL SPHERE  
============================================================

Sphere with same mass (36.4 kg):
  Radius: 7.97 cm

Running Monte Carlo simulation...
Criticality: 100%|████████████████████| Time: 0:01:48

Sphere k_eff = 1.321 ± 0.0040

Shape penalty (football vs sphere): 8.8%
```

| Geometry | k_eff | Reactivity |
|----------|-------|------------|
| Football | 1.214 | +17.6% |
| Sphere | 1.321 | +24.3% |

The football shape costs about **9% in reactivity** compared to an optimal sphere. The pointed ends cause extra neutron leakage.

## Critical Mass Analysis

How much plutonium is actually needed for criticality in each shape?

### Football Critical Mass

```
============================================================
CRITICAL MASS IN FOOTBALL SHAPE
============================================================

Searching for critical mass in range [5.0, 50.0] kg
Geometry: Football
Material: δ-phase Pu (Jezebel)

  Iteration  1: mass = 15.81 kg, k_eff = 0.9723 ± 0.0051
  Iteration  2: mass = 19.37 kg, k_eff = 1.0412 ± 0.0049
  Iteration  3: mass = 17.50 kg, k_eff = 1.0124 ± 0.0048
  Iteration  4: mass = 16.63 kg, k_eff = 0.9892 ± 0.0050
  Iteration  5: mass = 17.05 kg, k_eff = 1.0021 ± 0.0047

Critical mass (football shape): 17.1 kg
```

### Sphere Critical Mass

```
  Iteration  1: mass = 15.81 kg, k_eff = 1.0124 ± 0.0049
  Iteration  2: mass = 13.69 kg, k_eff = 0.9712 ± 0.0051
  Iteration  3: mass = 14.72 kg, k_eff = 0.9923 ± 0.0050
  Iteration  4: mass = 15.25 kg, k_eff = 1.0021 ± 0.0048

Critical mass (sphere): 15.3 kg
```

### Summary

| Shape | Critical Mass | Excess in Football |
|-------|--------------|-------------------|
| Sphere | 15.3 kg | — |
| Football | 17.1 kg | +12% |

The football requires **12% more material** to reach criticality. But at 36.4 kg, it contains **more than twice** the critical mass.

## Why Is the Football Supercritical?

Several factors contribute:

### 1. Abundant Fissile Material
The 36.4 kg of plutonium is well above critical mass for any geometry. Even the inefficient football shape can't waste enough neutrons to stay subcritical.

### 2. Favorable Density
δ-phase plutonium at 15.8 g/cm³ is dense enough that neutrons undergo multiple collisions before escaping. The mean free path (~4 cm) is much smaller than the football dimensions (~28 cm).

### 3. High k_∞ for Pu-239
Pure Pu-239 has k_∞ ≈ 2.9—every neutron absorbed produces nearly 3 new ones (on average). Even with significant leakage, k_eff > 1 is achievable.

## Physical Implications

### What Would Happen?

If you somehow assembled a Pu-filled football, several outcomes are possible:

**Scenario 1: Slow Assembly**
If assembled gradually, the system would go critical before completion. Alpha particles from decay and spontaneous fission neutrons would trigger the chain reaction. The result: a partial nuclear excursion that disperses the material (a "fizzle").

**Scenario 2: Fast Assembly (Impossible)**
Even at implosion speeds (~km/s), the high spontaneous fission rate of Pu-240 (~400,000 n/s per kg of WGPu) means predetonation is certain. The system would experience a low-yield explosion before reaching optimal configuration.

**Scenario 3: Hypothetical Instantaneous Assembly**
If magically assembled instantaneously, the 80-dollar supercritical system would produce a nuclear yield. The energy release:

$$E \approx \frac{1}{2}m_{fissioned}c^2 \cdot (efficiency)$$

With perhaps 1-5% of the material fissioning before disassembly, the yield would be in the kiloton range.

### Radiation Hazards

Even subcritical, a 36 kg Pu mass would be intensely radioactive:
- Alpha particles (short range, dangerous if inhaled)
- Gamma rays from decay products
- Neutrons from spontaneous fission

Handling would require specialized facilities (glove boxes, shielding, criticality safety protocols).

## Conclusions

| Question | Answer |
|----------|--------|
| Would a Pu-filled NFL football be critical? | **YES** |
| By how much? | k = 1.21, ρ = +18%, 80 dollars |
| How does shape affect it? | Football requires 12% more mass than sphere |
| Mass in football? | 36.4 kg (~2.2× critical mass) |
| Time scale? | Doubling every 47 ns |

The nuclear football is not just critical—it's deeply supercritical. The pointed shape causes extra leakage, but there's simply too much plutonium for that to matter.

## The Final Verdict

```
╔══════════════════════════════════════════════════════════╗
║                     FINAL VERDICT                        ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║   ✓ YES, a plutonium-filled NFL football would be a     ║
║     supercritical nuclear device.                        ║
║                                                          ║
║   The pointed shape is suboptimal but the ~36 kg of     ║
║   fissile material more than compensates.                ║
║                                                          ║
║   k_eff ≈ 1.21 (prompt supercritical)                   ║
║   Critical mass in football shape: ~17 kg               ║
║   Actual mass: 36 kg (2.1× critical)                    ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝
```

Don't try this at home.
