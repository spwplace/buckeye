# The Nucleus and Nuclear Forces

At the heart of every atom lies the nucleus—a tiny, dense core containing protons and neutrons bound together by the strong nuclear force. Understanding the nucleus is essential for understanding why some atoms split and release energy.

## The Scale of Things

The atom is mostly empty space. If the nucleus were a marble at the center of a football field, the electrons would be gnats buzzing around the parking lot. More precisely:

- **Atomic radius**: ~10⁻¹⁰ m (1 Ångström)
- **Nuclear radius**: ~10⁻¹⁵ m (1 femtometer)

The nucleus is 100,000 times smaller than the atom, yet contains 99.97% of its mass.

The nuclear radius follows an empirical formula:

$$R = r_0 A^{1/3}$$

where:
- \\(r_0 \approx 1.2\\) fm is a constant
- \\(A\\) is the mass number (protons + neutrons)

This cube-root scaling tells us something profound: nuclear density is constant. Double the nucleons, double the volume. The nucleus is like an incompressible liquid droplet.

## The Nuclear Force

What holds the nucleus together? Protons carry positive charge and should violently repel each other—the electrostatic force between two protons at nuclear distances is enormous:

$$F_\text{Coulomb} = \frac{e^2}{4\pi\epsilon_0 r^2} \approx 230 \text{ N at } r = 1 \text{ fm}$$

Yet nuclei exist. There must be an attractive force strong enough to overcome this repulsion.

This is the **strong nuclear force**, one of the four fundamental forces of nature. It has distinctive properties:

1. **Short range**: Essentially zero beyond ~2 fm
2. **Extremely strong**: About 100× stronger than electromagnetism at nuclear distances
3. **Charge-independent**: Acts equally between p-p, n-n, and p-n pairs
4. **Saturates**: Each nucleon only interacts with its nearest neighbors

The saturation property is key. Unlike gravity or electromagnetism (which are long-range), the strong force means each nucleon only "feels" the nucleons touching it. This is why nuclear density is constant—adding more nucleons doesn't compress the existing ones.

## Binding Energy

The mass of a nucleus is less than the sum of its constituent protons and neutrons. This "missing mass" has been converted to binding energy via Einstein's \\(E = mc^2\\).

The **binding energy** of a nucleus with \\(Z\\) protons and \\(N\\) neutrons (\\(A = Z + N\\)) is:

$$B(Z, N) = \left[ Z m_p + N m_n - M(Z, N) \right] c^2$$

where:
- \\(m_p = 938.272\\) MeV/c² (proton mass)
- \\(m_n = 939.565\\) MeV/c² (neutron mass)
- \\(M(Z, N)\\) is the actual nuclear mass

The binding energy per nucleon, \\(B/A\\), tells us how tightly bound the nucleus is. It varies across the periodic table:

| Element | A | B/A (MeV) |
|---------|---|-----------|
| ²H (deuterium) | 2 | 1.11 |
| ⁴He (helium) | 4 | 7.07 |
| ¹²C (carbon) | 12 | 7.68 |
| ⁵⁶Fe (iron) | 56 | 8.79 |
| ²³⁵U (uranium) | 235 | 7.59 |
| ²³⁹Pu (plutonium) | 239 | 7.56 |

The curve peaks around iron-56 at 8.79 MeV per nucleon. This has profound implications:

- **Fusion**: Light nuclei can release energy by combining (moving up the curve)
- **Fission**: Heavy nuclei can release energy by splitting (moving up the curve)

Plutonium-239 releases about 200 MeV per fission event, or roughly 1 MeV per nucleon.

## The Semi-Empirical Mass Formula

The liquid drop model provides a semi-empirical formula for nuclear binding energy:

$$B(Z, N) = a_V A - a_S A^{2/3} - a_C \frac{Z(Z-1)}{A^{1/3}} - a_A \frac{(N-Z)^2}{A} + \delta(A, Z)$$

Each term has a physical interpretation:

### Volume Term: \\(a_V A\\)
Binding is proportional to the number of nucleons (each contributes equally due to saturation). \\(a_V \approx 15.8\\) MeV.

### Surface Term: \\(-a_S A^{2/3}\\)
Nucleons at the surface have fewer neighbors and are less bound. Surface area scales as \\(A^{2/3}\\). \\(a_S \approx 18.3\\) MeV.

### Coulomb Term: \\(-a_C \frac{Z(Z-1)}{A^{1/3}}\\)
Protons repel each other. The total Coulomb energy of a uniformly charged sphere scales as \\(Z^2/R \propto Z^2/A^{1/3}\\). \\(a_C \approx 0.71\\) MeV.

### Asymmetry Term: \\(-a_A \frac{(N-Z)^2}{A}\\)
Nuclei prefer \\(N \approx Z\\) due to the Pauli exclusion principle. Deviation costs energy. \\(a_A \approx 23.2\\) MeV.

### Pairing Term: \\(\delta(A, Z)\\)
Nucleons prefer to pair up. Even-even nuclei are more stable:
$$\delta = \begin{cases} +a_P A^{-1/2} & \text{even-even} \\ 0 & \text{odd-A} \\ -a_P A^{-1/2} & \text{odd-odd} \end{cases}$$
with \\(a_P \approx 12\\) MeV.

## The Valley of Stability

Plotting stable nuclei on an \\(N\\) vs \\(Z\\) chart reveals the **valley of stability**—a band where nuclei are stable against radioactive decay. Key features:

1. **Light nuclei**: \\(N \approx Z\\) (the asymmetry term dominates)
2. **Heavy nuclei**: \\(N > Z\\) (extra neutrons dilute Coulomb repulsion)
3. **Beyond lead**: No stable nuclei exist for \\(Z > 82\\)
4. **Drip lines**: Limits where nuclei become unbound

Plutonium-239 (\\(Z = 94\\), \\(N = 145\\)) is well beyond the valley of stability. It's radioactive, decaying via alpha emission with a half-life of 24,100 years. But for our purposes, it's stable enough—and critically, it's fissile.

## Nuclear Instability

Heavy nuclei are unstable for two reasons:

1. **Coulomb repulsion** grows as \\(Z^2\\) while the strong force only grows as \\(A\\)
2. **Surface energy** of heavy nuclei makes them easier to deform

The liquid drop model predicts a **fissility parameter**:

$$x = \frac{E_\text{Coulomb}}{2 E_\text{Surface}} = \frac{Z^2/A}{50.88}$$

When \\(x \geq 1\\), the Coulomb repulsion exceeds the surface tension and the nucleus spontaneously fissions. For ²³⁹Pu:

$$x = \frac{94^2/239}{50.88} = 0.73$$

Plutonium is below the spontaneous fission threshold but close enough that a small energy input (from an absorbed neutron) can push it over the edge.

## Key Takeaways

1. Nuclei are held together by the strong nuclear force, which is short-range and saturates
2. Binding energy per nucleon peaks at iron-56, making both fusion (light → medium) and fission (heavy → medium) energetically favorable
3. Heavy nuclei are inherently unstable due to Coulomb repulsion
4. Plutonium-239 is close to the spontaneous fission threshold, making it easily fissile

In the next chapter, we'll see exactly what happens when a neutron strikes a plutonium nucleus.
