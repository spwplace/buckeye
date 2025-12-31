# Neutron Interactions

Before a neutron can cause fission, it must interact with a nucleus. Understanding these interactions is essential for predicting how neutrons travel through matter.

## Types of Neutron Reactions

When a neutron approaches a nucleus, several things can happen:

### Elastic Scattering (n, n)
The neutron bounces off the nucleus like a billiard ball. Both kinetic energy and momentum are conserved. The nucleus recoils but remains in its ground state.

$$n + A \rightarrow n + A$$

This is the most common interaction for fast neutrons in heavy materials.

### Inelastic Scattering (n, n')
The neutron excites the nucleus to a higher energy state. The nucleus later de-excites by emitting gamma rays. The neutron loses more energy than in elastic scattering.

$$n + A \rightarrow n' + A^* \rightarrow n' + A + \gamma$$

This requires the neutron to have enough energy to excite the nucleus (threshold ~40 keV for heavy nuclei).

### Radiative Capture (n, γ)
The neutron is absorbed, and the binding energy is released as gamma rays:

$$n + A \rightarrow (A+1) + \gamma$$

The product nucleus may be radioactive. For Pu-239:
$$n + {}^{239}\text{Pu} \rightarrow {}^{240}\text{Pu} + \gamma$$

Capture removes neutrons from the chain reaction.

### Fission (n, f)
Discussed in the previous chapter. The nucleus absorbs the neutron and splits:

$$n + {}^{239}\text{Pu} \rightarrow \text{fragments} + 2\text{-}3 n + \gamma + 200 \text{ MeV}$$

This produces neutrons—the engine of the chain reaction.

### Other Reactions
Less common for our purposes:
- **(n, 2n)**: Neutron knocks out two neutrons
- **(n, p)**: Neutron knocks out a proton
- **(n, α)**: Neutron causes alpha emission

## The Cross Section Concept

The probability of a reaction is quantified by its **cross section**, denoted σ. Think of it as the effective "target area" the nucleus presents to the neutron.

If a nucleus had geometric area \\(\pi R^2\\), and every neutron hitting it reacted, that would be the cross section. In reality, quantum mechanics makes cross sections very different from geometric areas—they can be much larger (resonances) or much smaller (threshold reactions).

### Units: The Barn

Cross sections are measured in **barns**:
$$1 \text{ barn} = 10^{-24} \text{ cm}^2$$

The name comes from Manhattan Project physicists who thought nuclear cross sections were "as big as a barn" compared to the geometric nuclear area (~10⁻²⁵ cm²).

Typical values:
- Pu-239 fission: ~1.7 barns (1 MeV neutrons)
- Pu-239 elastic scattering: ~4 barns (1 MeV)
- Pu-239 capture: ~0.01 barns (1 MeV)

### Microscopic vs Macroscopic

The **microscopic cross section** σ refers to a single nucleus. The **macroscopic cross section** Σ accounts for the number density of nuclei:

$$\Sigma = N \sigma$$

where N is the number density (atoms/cm³). Macroscopic cross sections have units of cm⁻¹ and represent the probability of interaction per unit path length.

For Pu-239 at density 15.8 g/cm³:
$$N = \frac{\rho N_A}{A} = \frac{15.8 \times 6.02 \times 10^{23}}{239} = 3.98 \times 10^{22} \text{ atoms/cm}^3$$

The macroscopic total cross section at 1 MeV:
$$\Sigma_t = (3.98 \times 10^{22})(6.2 \times 10^{-24}) = 0.247 \text{ cm}^{-1}$$

This means a neutron travels about 4 cm on average before interacting.

## Mean Free Path

The **mean free path** is the average distance a neutron travels between collisions:

$$\lambda = \frac{1}{\Sigma_t}$$

For the example above, λ ≈ 4 cm. In a plutonium football (longest dimension ~28 cm), a neutron will undergo roughly 7 collisions before escaping.

## Reaction Rates

The rate at which neutrons undergo a particular reaction is:

$$R = \Sigma \phi$$

where φ is the **neutron flux** (neutrons per cm² per second). For a single neutron with velocity v:

$$\phi = nv$$

where n is the neutron density. The reaction rate can also be written:

$$R = N \sigma v n = N \langle \sigma v \rangle n$$

where ⟨σv⟩ is the thermally-averaged reaction rate.

## Energy Dependence

Cross sections vary dramatically with neutron energy. This is why we need continuous-energy Monte Carlo rather than simplified "one-group" calculations.

### 1/v Region (Thermal)

At low energies (E < 1 eV), many cross sections follow:

$$\sigma(E) \propto \frac{1}{v} \propto \frac{1}{\sqrt{E}}$$

This is because slow neutrons spend more time near the nucleus, giving the reaction more time to occur.

### Resonance Region (Epithermal)

Between ~1 eV and ~10 keV, cross sections exhibit sharp peaks called **resonances**. These occur when the neutron's energy matches an excited state of the compound nucleus.

The resonance shape is given by the **Breit-Wigner formula**:

$$\sigma(E) = \sigma_0 \frac{\Gamma^2/4}{(E - E_0)^2 + \Gamma^2/4}$$

where:
- E₀ is the resonance energy
- Γ is the resonance width
- σ₀ is the peak cross section

Heavy nuclei like Pu-239 have thousands of resonances in the epithermal region. Resolving these accurately is critical for reactor calculations.

### Fast Region

Above ~10 keV, resonances begin to overlap and the cross section becomes smoother. In this region, optical model calculations provide accurate predictions.

Key features for Pu-239:
- Fission cross section relatively flat at ~1.6-2.0 barns
- Elastic scattering decreases slowly from ~5 to ~3 barns
- Inelastic scattering rises from threshold, peaks around 2 MeV
- Capture falls rapidly from the thermal value

## Scattering Kinematics

When a neutron scatters elastically off a nucleus, it loses energy. The kinematics are determined by conservation of momentum and energy.

### Energy Transfer

In the center-of-mass frame, elastic scattering can be described by the scattering angle μ_cm = cos(θ_cm). The relationship between initial and final neutron energies in the lab frame is:

$$E' = E \cdot \frac{A^2 + 2A\mu_{cm} + 1}{(A+1)^2}$$

where A is the target mass number.

The minimum energy after collision (backscatter, μ_cm = -1):

$$E'_\min = \alpha E$$

where:

$$\alpha = \left(\frac{A-1}{A+1}\right)^2$$

For Pu-239 (A = 239):
$$\alpha = \left(\frac{238}{240}\right)^2 = 0.983$$

A neutron loses at most 1.7% of its energy per collision with plutonium. This is why heavy nuclei are poor moderators—it takes many collisions to slow neutrons down.

Compare to hydrogen (A = 1):
$$\alpha = 0$$

A neutron can lose all its energy in a single collision with hydrogen. This is why water is an effective moderator.

### Average Energy Loss

The average logarithmic energy loss per collision is:

$$\xi = 1 + \frac{\alpha \ln\alpha}{1-\alpha}$$

For heavy nuclei, this simplifies to:
$$\xi \approx \frac{2}{A + 2/3}$$

For Pu-239: ξ ≈ 0.0084, meaning each collision reduces ln(E) by about 0.8%.

### Angular Distribution

The angular distribution in the center-of-mass frame is often nearly isotropic for heavy nuclei:

$$P(\mu_{cm}) = \frac{1}{2}, \quad -1 \leq \mu_{cm} \leq 1$$

The transformation to lab frame:

$$\mu_{lab} = \frac{1 + A\mu_{cm}}{\sqrt{A^2 + 2A\mu_{cm} + 1}}$$

For heavy nuclei, μ_lab ≈ μ_cm ≈ 0 on average (isotropic scattering in both frames).

## Cross Section Data

Real calculations require tabulated cross section data. The standard source is the **Evaluated Nuclear Data File (ENDF)**:

- **ENDF/B-VIII.0**: US evaluation (2018)
- **JEFF-3.3**: European evaluation (2017)
- **JENDL-5**: Japanese evaluation (2021)

These contain:
- Point-wise cross sections on fine energy grids
- Resonance parameters for exact reconstruction
- Angular distributions
- Fission spectrum parameters
- Delayed neutron data

Our code uses simplified tabulations derived from ENDF/B-VIII.0, interpolated via log-log interpolation (standard practice for cross sections).

## Summary

| Reaction | Effect on Neutrons | Cross Section |
|----------|-------------------|---------------|
| Elastic scatter | Changes direction, small E loss | ~4 barns |
| Inelastic scatter | Changes direction, larger E loss | ~2 barns |
| Fission | Removes 1, creates 2-3 | ~1.7 barns |
| Capture | Removes 1 | ~0.01 barns |

The competition between fission (creates neutrons) and absorption+leakage (removes neutrons) determines whether a chain reaction can sustain itself—the subject of the next chapter.
