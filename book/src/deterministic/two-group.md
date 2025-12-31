# Two-Group Theory

One-group diffusion theory averages all neutron energies together. This loses important physics in systems where neutrons span a wide energy range. **Two-group theory** provides a middle ground between one-group and continuous-energy treatment.

## Energy Group Structure

The neutron energy spectrum spans ~10 decades (10⁻⁵ to 10⁷ eV). We divide this into two groups:

| Group | Energy Range | Character | Typical Energy |
|-------|-------------|-----------|----------------|
| Group 1 (Fast) | E > 1 MeV | Fission spectrum | 2 MeV |
| Group 2 (Epithermal) | 0.1 eV - 1 MeV | Slowed neutrons | 0.4 MeV |

For fast systems like bare plutonium, most neutrons remain in Group 1, but down-scattering to Group 2 affects the overall neutron balance.

## Two-Group Equations

The steady-state two-group diffusion equations are:

**Group 1 (Fast):**
$$-D_1\nabla^2\phi_1 + \Sigma_{r1}\phi_1 = \chi_1(\nu_1\Sigma_{f1}\phi_1 + \nu_2\Sigma_{f2}\phi_2)$$

**Group 2 (Epithermal):**
$$-D_2\nabla^2\phi_2 + \Sigma_{a2}\phi_2 = \chi_2(\nu_1\Sigma_{f1}\phi_1 + \nu_2\Sigma_{f2}\phi_2) + \Sigma_{s,1\to2}\phi_1$$

where:
- Σ_r1 = Σ_a1 + Σ_s,1→2 (removal from Group 1)
- χ_1, χ_2 = fission spectrum fractions (χ_1 ≈ 0.75, χ_2 ≈ 0.25)
- Σ_s,1→2 = down-scatter cross section

## Two-Group Cross Sections

The Python code defines group constants for Pu-239:

```python
PU239_TWO_GROUP = TwoGroupData(
    name="Pu-239",
    A=239.0521634,
    rho_0=19.86e3,
    group1=GroupConstants(
        sigma_f=1.96,      # fission (barns)
        sigma_c=0.045,     # capture
        sigma_s=5.4,       # scattering
        sigma_tr=5.6,      # transport
        nu=3.12,           # neutrons/fission
        E_avg=2.0,         # average energy (MeV)
    ),
    group2=GroupConstants(
        sigma_f=1.78,
        sigma_c=0.16,
        sigma_s=6.2,
        sigma_tr=6.5,
        nu=2.94,
        E_avg=0.4,
    ),
    sigma_s12=0.12,        # down-scatter
)
```

Note that:
- Fission cross section is similar in both groups
- Capture is much higher in Group 2 (1/v behavior)
- Down-scatter cross section is small for heavy nuclei

## k_∞ Calculation

For an infinite medium, the two-group k_∞ requires solving an eigenvalue problem:

$$\mathbf{L}\boldsymbol{\phi} = \frac{1}{k}\mathbf{F}\boldsymbol{\phi}$$

where:
$$\mathbf{L} = \begin{pmatrix} \Sigma_{r1} & 0 \\ -\Sigma_{s12} & \Sigma_{a2} \end{pmatrix}, \quad
\mathbf{F} = \begin{pmatrix} \chi_1\nu_1\Sigma_{f1} & \chi_1\nu_2\Sigma_{f2} \\ \chi_2\nu_1\Sigma_{f1} & \chi_2\nu_2\Sigma_{f2} \end{pmatrix}$$

The dominant eigenvalue of L⁻¹F gives k_∞.

```python
def compute_two_group_k_inf(xs: TwoGroupMacroXS) -> Tuple[float, Dict]:
    Sigma_r1 = xs.Sigma_a1 + xs.Sigma_s12
    Sigma_r2 = xs.Sigma_a2
    
    nu_Sigma_f1 = xs.nu1 * xs.Sigma_f1
    nu_Sigma_f2 = xs.nu2 * xs.Sigma_f2
    
    # Build L⁻¹F matrix
    L_inv_11 = 1.0 / Sigma_r1
    L_inv_21 = xs.Sigma_s12 / (Sigma_r1 * Sigma_r2)
    L_inv_22 = 1.0 / Sigma_r2
    
    A_11 = L_inv_11 * CHI_1 * nu_Sigma_f1
    A_12 = L_inv_11 * CHI_1 * nu_Sigma_f2
    A_21 = L_inv_21 * CHI_1 * nu_Sigma_f1 + L_inv_22 * CHI_2 * nu_Sigma_f1
    A_22 = L_inv_21 * CHI_1 * nu_Sigma_f2 + L_inv_22 * CHI_2 * nu_Sigma_f2
    
    # Eigenvalues of 2×2 matrix
    trace = A_11 + A_22
    det = A_11 * A_22 - A_12 * A_21
    k_inf = (trace + np.sqrt(trace**2 - 4*det)) / 2
    
    return k_inf, {...}
```

## Two-Group Critical Radius

With geometric buckling B²:

$$\mathbf{L}_B\boldsymbol{\phi} = \frac{1}{k}\mathbf{F}\boldsymbol{\phi}$$

where L_B includes leakage terms:
$$\mathbf{L}_B = \begin{pmatrix} \Sigma_{r1} + D_1B^2 & 0 \\ -\Sigma_{s12} & \Sigma_{a2} + D_2B^2 \end{pmatrix}$$

Criticality occurs when k = 1, which determines B² and hence R_c.

## Time-Dependent Two-Group

For kinetics, we track both group populations:

$$\frac{dn_1}{dt} = \chi_1 S - (\Sigma_{a1} + \Sigma_{s12})v_1 n_1$$
$$\frac{dn_2}{dt} = \chi_2 S + \Sigma_{s12}v_1 n_1 - \Sigma_{a2}v_2 n_2$$

where S is the total fission source.

This gives a 2×2 eigenvalue problem for the system alpha:

$$\frac{d\mathbf{n}}{dt} = \mathbf{A}\mathbf{n}$$

The dominant eigenvalue of **A** is the effective α.

## Group Coupling Effects

In fast systems, down-scattering is weak (Σ_s12 ≈ 0.12 barns), so:
- Most neutrons remain in Group 1
- Group 2 population is ~10-20% of Group 1
- k_∞ is dominated by Group 1 fission

The two-group model captures:
1. Spectral hardening during compression
2. Different fission probabilities at different energies
3. Energy-dependent leakage

## Comparison: One-Group vs Two-Group

| Property | One-Group | Two-Group |
|----------|-----------|-----------|
| k_∞ (Pu-239) | 2.82 | 2.74 |
| Critical radius | 4.9 cm | 5.0 cm |
| Critical mass | 10.0 kg | 10.4 kg |

The two-group model gives ~4% higher critical mass due to better treatment of down-scatter losses.

## Why Not More Groups?

Industrial reactor codes use 2-70 energy groups. Why stop at two?

For fast systems:
- Spectrum is relatively flat above 100 keV
- No thermal peak to resolve
- Two groups capture essential physics
- More groups add computation without proportional accuracy gain

For thermal reactors, many groups are needed to:
- Resolve resonances
- Track thermalization
- Capture poison effects

## Python Implementation

The two-group simulation (`two_group_physics.py`) includes:

```python
class TwoGroupSimulation:
    def simulate(self, n1_0, n2_0, t_max):
        # Initial flux ratio from eigenvalue problem
        xs_0 = compute_two_group_xs(self.core_data, self.compression)
        _, eigvecs, _ = compute_two_group_alpha(xs_0, self.R_0)
        
        # Time-stepping with BDF method
        sol = solve_ivp(rhs, (0, t_max), y0, method="BDF", ...)
        
        return results
```

The simulation tracks both group densities and couples to hydrodynamics.

## Results

For 6.2 kg Pu-239 at 2.5× compression:

```
Initial conditions (two-group):
  R_0 = 4.17 cm
  compression = 2.50
  α_dominant = 2.1×10⁸ /s
  n1_0/n2_0 = 1.0×10¹⁰ / 2.1×10⁹

Final state:
  t = 0.38 μs
  Yield = 11.8 kt
```

The two-group yield is ~5% lower than one-group due to better loss treatment.

## Summary

Two-group theory provides:
- Energy-dependent physics without continuous-energy complexity
- Improved accuracy over one-group (~5-10%)
- Physical insight into spectral effects
- Manageable analytical framework

For fast bare systems, two groups capture the essential physics. For more complex problems (moderated systems, heterogeneous cores), more energy groups may be needed.
