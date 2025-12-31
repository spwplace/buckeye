# Nuclear Fission

When a heavy nucleus like plutonium-239 absorbs a neutron, it becomes unstable and splits into two lighter nuclei, releasing enormous energy and additional neutrons. This is nuclear fission—the process that powers reactors and bombs.

## The Fission Process

Fission occurs in stages:

### 1. Neutron Absorption
A neutron collides with a ²³⁹Pu nucleus and is captured, forming a **compound nucleus**:

$$n + {}^{239}\text{Pu} \rightarrow {}^{240}\text{Pu}^*$$

The asterisk indicates an excited state. The compound nucleus has absorbed the neutron's kinetic energy plus its binding energy in the new configuration (about 6.5 MeV for Pu-239).

### 2. Nuclear Deformation
The compound nucleus begins oscillating like a liquid droplet. The excitation energy causes large-amplitude vibrations, stretching the nucleus into an elongated shape.

The competition between surface tension (trying to restore spherical shape) and Coulomb repulsion (trying to push the protons apart) determines the nucleus's fate.

### 3. Scission
If the deformation is large enough, the nucleus reaches a **saddle point**—a point of no return where Coulomb repulsion wins. The nucleus necks down and splits into two fragments:

$${}^{240}\text{Pu}^* \rightarrow {}^{A_1}\text{X} + {}^{A_2}\text{Y} + \nu n + \gamma$$

where:
- X and Y are fission fragments (typically unequal mass)
- $\nu$ neutrons are released (typically 2-3)
- Gamma rays carry away additional energy

### 4. Fragment Acceleration
The fission fragments are born close together with enormous Coulomb repulsion. They accelerate apart, converting ~170 MeV of Coulomb potential energy into kinetic energy.

## Energy Release

Fission releases roughly 200 MeV per event—about a million times more than chemical reactions. The energy breakdown:

| Component | Energy (MeV) |
|-----------|-------------|
| Fission fragment kinetic energy | ~165 |
| Prompt neutron kinetic energy | ~5 |
| Prompt gamma rays | ~7 |
| Fission fragment beta decay | ~8 |
| Fission fragment gamma decay | ~7 |
| Neutrinos (lost) | ~12 |
| **Total** | **~204** |

The "useful" energy (excluding neutrinos, which escape) is about 192 MeV per fission. For plutonium-239, this corresponds to:

$$E = 192 \text{ MeV/fission} \times \frac{6.02 \times 10^{23} \text{ atoms/mol}}{239 \text{ g/mol}} = 7.7 \times 10^{13} \text{ J/kg}$$

Complete fission of 1 kg of Pu-239 releases the energy equivalent of about 18,000 tons of TNT.

## Fission Products

Fission produces a characteristic distribution of fragment masses. For thermal neutron fission of Pu-239:

- **Light fragment peak**: A ≈ 100 (Mo, Tc, Ru)
- **Heavy fragment peak**: A ≈ 138 (Ba, La, Ce)
- **Symmetric fission**: Rare (valley between peaks)

The asymmetry arises from nuclear shell effects—fragments near "magic numbers" (filled nuclear shells) are energetically favored.

Fission products are typically neutron-rich and undergo beta decay chains until reaching stability. Some important fission products:

- **¹³⁷Cs** (cesium-137): 30-year half-life, major radioactive waste concern
- **¹³¹I** (iodine-131): 8-day half-life, thyroid hazard
- **⁹⁰Sr** (strontium-90): 29-year half-life, bone-seeking
- **¹³⁵Xe** (xenon-135): Neutron poison in reactors

## Prompt Neutrons

The neutrons released during fission are crucial for chain reactions. Key parameters:

### Average Number: ν̄

The average number of neutrons per fission depends on the fissioning nucleus and the incident neutron energy:

$$\bar{\nu}(E) = \bar{\nu}_0 + \alpha E$$

For Pu-239:
- $\bar{\nu}_0 = 2.874$ (at thermal energies)
- $\alpha = 0.148$ MeV⁻¹

At 1 MeV, $\bar{\nu} \approx 3.02$ neutrons per fission.

### Energy Spectrum: χ(E)

Fission neutrons are born with a characteristic energy distribution called the **fission spectrum** or **Watt spectrum**:

$$\chi(E) \propto e^{-E/a} \sinh\sqrt{bE}$$

For Pu-239:
- $a = 0.966$ MeV
- $b = 2.842$ MeV⁻¹

This gives:
- **Most probable energy**: ~0.7 MeV
- **Average energy**: ~2.0 MeV
- **High-energy tail**: Extends beyond 10 MeV

### Angular Distribution

Prompt fission neutrons are emitted **isotropically** in the center-of-mass frame (which is essentially the lab frame for heavy nuclei).

## Delayed Neutrons

A small fraction (~0.2-0.7%) of fission neutrons are **delayed**—emitted seconds to minutes after fission from excited fission product nuclei. These are crucial for reactor control:

| Group | Half-life | Yield (Pu-239) |
|-------|-----------|----------------|
| 1 | 54.5 s | 0.000086 |
| 2 | 21.8 s | 0.000637 |
| 3 | 5.98 s | 0.000491 |
| 4 | 2.23 s | 0.000743 |
| 5 | 0.50 s | 0.000240 |
| 6 | 0.18 s | 0.000087 |

The total delayed neutron fraction for Pu-239 is β ≈ 0.0022 (0.22%). Compare to U-235 with β ≈ 0.0065.

This matters because:
- **Prompt critical**: k = 1 using only prompt neutrons (instant exponential growth)
- **Delayed critical**: k = 1 including delayed neutrons (controllable growth)

The small β for plutonium makes Pu reactors harder to control than uranium reactors.

## Spontaneous Fission

Even without neutron bombardment, heavy nuclei can spontaneously tunnel through the fission barrier. The spontaneous fission rate varies dramatically:

| Nuclide | SF half-life |
|---------|-------------|
| ²³⁵U | 3.5 × 10¹⁷ years |
| ²³⁸U | 8.2 × 10¹⁵ years |
| ²³⁹Pu | 8.0 × 10¹⁵ years |
| ²⁴⁰Pu | 1.14 × 10¹¹ years |

Note that Pu-240 has a spontaneous fission rate 70,000 times higher than Pu-239. This is why weapons-grade plutonium must have low Pu-240 content—too much causes predetonation.

A bare 10 kg sphere of weapons-grade plutonium (6% Pu-240) experiences about 400,000 spontaneous fission events per second. Each provides a neutron that could start a chain reaction if the assembly is supercritical.

## Fissile vs Fissionable

Important terminology:

- **Fissile**: Can sustain a chain reaction with neutrons of any energy. Examples: U-233, U-235, Pu-239, Pu-241.

- **Fissionable**: Can fission, but only with high-energy neutrons (above threshold). Examples: U-238, Th-232.

- **Fertile**: Can be converted to fissile material by neutron capture. Examples: U-238 → Pu-239, Th-232 → U-233.

Our nuclear football is filled with Pu-239, which is fissile—it will sustain a chain reaction regardless of neutron energy.

## The Chain Reaction

Combine these facts:
1. One fission releases ~3 neutrons
2. Each neutron can cause another fission
3. Each fission releases enormous energy

If, on average, more than one neutron from each fission causes another fission, the neutron population grows exponentially. This is a **chain reaction**.

The parameter controlling this is **k**, the multiplication factor—the subject of our next chapter.

## Key Takeaways

1. Fission splits a heavy nucleus into two fragments, releasing ~200 MeV and 2-3 neutrons
2. The neutron energy spectrum follows the Watt distribution, peaking around 0.7 MeV
3. Plutonium-239 is fissile (sustains chain reactions at any neutron energy)
4. Delayed neutrons are crucial for reactor control but Pu-239 has few of them
5. Pu-240 impurity causes spontaneous fission, potentially triggering predetonation
