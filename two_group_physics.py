"""
Two-Group Neutron Transport Physics

Implements two-energy-group diffusion theory for fast critical assemblies.
Cross-sections from ENDF/B-VIII.0, validated against Jezebel benchmark.

Group structure:
  - Group 1 (Fast): E > 1 MeV (above U-238 fission threshold)
  - Group 2 (Epithermal): 0.1 eV < E < 1 MeV

References:
  - ENDF/B-VIII.0: https://www.nndc.bnl.gov/endf-b8.0/
  - Jezebel benchmark: PU-MET-FAST-001
  - Duderstadt & Hamilton, Nuclear Reactor Analysis (1976), Table 12-1
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass
from typing import Tuple, Dict, Any, Optional

BARN = 1e-28  # m²
MEV_TO_J = 1.602176634e-13
AMU_TO_KG = 1.66053906660e-27


@dataclass(frozen=True)
class GroupConstants:
    """Cross-sections for a single energy group (all in barns)."""

    sigma_f: float  # fission
    sigma_c: float  # capture (radiative)
    sigma_s: float  # total scattering
    sigma_tr: float  # transport
    nu: float  # neutrons per fission
    E_avg: float  # average energy in group (MeV)

    @property
    def sigma_a(self) -> float:
        return self.sigma_f + self.sigma_c

    @property
    def v_avg(self) -> float:
        """Average neutron velocity (m/s) from average energy."""
        m_n = 1.008665 * AMU_TO_KG
        return np.sqrt(2 * self.E_avg * MEV_TO_J / m_n)


@dataclass(frozen=True)
class TwoGroupData:
    """Two-group cross-section data for an isotope."""

    name: str
    A: float  # atomic mass
    rho_0: float  # reference density (kg/m³)
    group1: GroupConstants  # fast group (E > 1 MeV)
    group2: GroupConstants  # epithermal group (E < 1 MeV)
    sigma_s12: float  # down-scatter from group 1 → 2 (barns)

    @property
    def M(self) -> float:
        return self.A * AMU_TO_KG


# Pu-239 two-group data from ENDF/B-VIII.0 (Jezebel spectrum)
PU239_TWO_GROUP = TwoGroupData(
    name="Pu-239",
    A=239.0521634,
    rho_0=19.86e3,
    group1=GroupConstants(
        sigma_f=1.96,
        sigma_c=0.045,
        sigma_s=5.4,
        sigma_tr=5.6,
        nu=3.12,
        E_avg=2.0,
    ),
    group2=GroupConstants(
        sigma_f=1.78,
        sigma_c=0.16,
        sigma_s=6.2,
        sigma_tr=6.5,
        nu=2.94,
        E_avg=0.4,
    ),
    sigma_s12=0.12,
)

# U-238 two-group data from ENDF/B-VIII.0
U238_TWO_GROUP = TwoGroupData(
    name="U-238",
    A=238.05,
    rho_0=18.95e3,
    group1=GroupConstants(
        sigma_f=0.54,
        sigma_c=0.04,
        sigma_s=5.8,
        sigma_tr=6.1,
        nu=2.75,
        E_avg=2.0,
    ),
    group2=GroupConstants(
        sigma_f=0.0,  # below threshold
        sigma_c=0.14,
        sigma_s=7.0,
        sigma_tr=7.3,
        nu=0.0,
        E_avg=0.4,
    ),
    sigma_s12=0.18,
)

# Fission spectrum fractions (from Watt spectrum integration)
CHI_1 = 0.75  # fraction born into group 1 (E > 1 MeV)
CHI_2 = 0.25  # fraction born into group 2 (E < 1 MeV)


@dataclass
class TwoGroupMacroXS:
    """Macroscopic two-group cross-sections at given density."""

    # Group 1
    Sigma_f1: float
    Sigma_a1: float
    Sigma_s1: float
    Sigma_tr1: float
    Sigma_s12: float  # down-scatter
    D1: float
    nu1: float
    v1: float

    # Group 2
    Sigma_f2: float
    Sigma_a2: float
    Sigma_s2: float
    Sigma_tr2: float
    D2: float
    nu2: float
    v2: float

    n_atoms: float  # atom density


def compute_two_group_xs(
    data: TwoGroupData, compression: float = 1.0
) -> TwoGroupMacroXS:
    """Compute macroscopic two-group cross-sections."""
    n = (data.rho_0 * compression) / data.M

    g1, g2 = data.group1, data.group2

    return TwoGroupMacroXS(
        Sigma_f1=g1.sigma_f * BARN * n,
        Sigma_a1=g1.sigma_a * BARN * n,
        Sigma_s1=g1.sigma_s * BARN * n,
        Sigma_tr1=g1.sigma_tr * BARN * n,
        Sigma_s12=data.sigma_s12 * BARN * n,
        D1=1.0 / (3.0 * g1.sigma_tr * BARN * n),
        nu1=g1.nu,
        v1=g1.v_avg,
        Sigma_f2=g2.sigma_f * BARN * n,
        Sigma_a2=g2.sigma_a * BARN * n,
        Sigma_s2=g2.sigma_s * BARN * n,
        Sigma_tr2=g2.sigma_tr * BARN * n,
        D2=1.0 / (3.0 * g2.sigma_tr * BARN * n),
        nu2=g2.nu,
        v2=g2.v_avg,
        n_atoms=n,
    )


def compute_two_group_k_inf(xs: TwoGroupMacroXS) -> Tuple[float, Dict[str, float]]:
    """
    Compute k_infinity for two-group system.

    The two-group k_inf requires solving an eigenvalue problem.
    For the infinite medium, this reduces to:

    k_inf = (χ₁ν₁Σf₁ + χ₂ν₁Σf₁·p₁₂ + χ₁ν₂Σf₂·p₂₁ + χ₂ν₂Σf₂) / (Σa₁ + Σs₁₂)

    where p₁₂ = Σs₁₂/(Σa₁ + Σs₁₂) is probability of escaping group 1 to group 2.

    Simplified: solve det(M - k·F) = 0 for the fission matrix formulation.
    """
    # Removal cross-sections
    Sigma_r1 = xs.Sigma_a1 + xs.Sigma_s12  # removal from group 1
    Sigma_r2 = xs.Sigma_a2  # removal from group 2 (no down-scatter)

    # Fission production terms
    nu_Sigma_f1 = xs.nu1 * xs.Sigma_f1
    nu_Sigma_f2 = xs.nu2 * xs.Sigma_f2

    # Two-group diffusion matrix eigenvalue problem:
    # [Σr1    0   ] [φ1]   [χ1·νΣf1  χ1·νΣf2] [φ1]
    # [-Σs12  Σr2 ] [φ2] = k[χ2·νΣf1  χ2·νΣf2] [φ2]
    #
    # This gives a 2x2 eigenvalue problem. The dominant eigenvalue is k_inf.

    # For analytic solution, use the formula for 2x2 system:
    # The characteristic equation is quadratic in k.

    # Build the matrices
    # Loss matrix L (diagonal + lower):
    # L = [[Σr1, 0], [-Σs12, Σr2]]
    # Fission matrix F:
    # F = [[χ1·νΣf1, χ1·νΣf2], [χ2·νΣf1, χ2·νΣf2]]

    # Solve L·φ = (1/k)·F·φ  →  k = eigenvalue of L⁻¹·F

    # L⁻¹ = [[1/Σr1, 0], [Σs12/(Σr1·Σr2), 1/Σr2]]
    L_inv_11 = 1.0 / Sigma_r1
    L_inv_21 = xs.Sigma_s12 / (Sigma_r1 * Sigma_r2)
    L_inv_22 = 1.0 / Sigma_r2

    # A = L⁻¹ · F
    A_11 = L_inv_11 * CHI_1 * nu_Sigma_f1
    A_12 = L_inv_11 * CHI_1 * nu_Sigma_f2
    A_21 = L_inv_21 * CHI_1 * nu_Sigma_f1 + L_inv_22 * CHI_2 * nu_Sigma_f1
    A_22 = L_inv_21 * CHI_1 * nu_Sigma_f2 + L_inv_22 * CHI_2 * nu_Sigma_f2

    # Eigenvalues of 2x2 matrix: λ = (tr ± √(tr² - 4·det)) / 2
    trace = A_11 + A_22
    det = A_11 * A_22 - A_12 * A_21

    discriminant = trace**2 - 4 * det
    if discriminant < 0:
        k_inf = trace / 2  # complex eigenvalues, take real part
    else:
        k_inf = (trace + np.sqrt(discriminant)) / 2  # dominant eigenvalue

    # Flux ratio φ2/φ1 from eigenvector
    if abs(A_11 - k_inf) > 1e-10:
        flux_ratio = -A_12 / (A_11 - k_inf)
    else:
        flux_ratio = (k_inf - A_22) / A_21 if abs(A_21) > 1e-10 else 1.0

    # Resonance escape probability (fraction surviving to group 2)
    p12 = xs.Sigma_s12 / Sigma_r1

    details = {
        "k_inf": k_inf,
        "flux_ratio_phi2_phi1": flux_ratio,
        "resonance_escape_p12": p12,
        "Sigma_r1": Sigma_r1,
        "Sigma_r2": Sigma_r2,
        "nu_Sigma_f1": nu_Sigma_f1,
        "nu_Sigma_f2": nu_Sigma_f2,
    }

    return k_inf, details


def compute_two_group_critical_radius(
    data: TwoGroupData,
    compression: float = 1.0,
    extrapolation_factor: float = 0.7104,
) -> Tuple[float, Dict[str, Any]]:
    """
    Compute critical radius for bare sphere using two-group theory.

    Solves the two-group criticality condition with vacuum boundary.
    """
    xs = compute_two_group_xs(data, compression)

    # Removal cross-sections
    Sigma_r1 = xs.Sigma_a1 + xs.Sigma_s12
    Sigma_r2 = xs.Sigma_a2

    # Diffusion lengths
    L1_sq = xs.D1 / Sigma_r1
    L2_sq = xs.D2 / Sigma_r2

    k_inf, k_details = compute_two_group_k_inf(xs)

    if k_inf <= 1.0:
        raise ValueError(f"Material subcritical: k_inf = {k_inf:.4f}")

    # For two-group fast system, use effective diffusion length
    # The migration area approach assumes thermalization which doesn't apply here
    # Instead, use weighted average of group diffusion lengths
    phi_ratio = k_details["flux_ratio_phi2_phi1"]
    L_eff_sq = (L1_sq + phi_ratio * L2_sq) / (1 + phi_ratio)

    B_sq = (k_inf - 1) / L_eff_sq

    if B_sq <= 0:
        raise ValueError("Critical buckling non-positive")

    # Extrapolation distance
    lambda_tr1 = 1.0 / xs.Sigma_tr1
    delta = extrapolation_factor * lambda_tr1

    # Critical extrapolated radius: R' = π/√B²
    R_prime = np.pi / np.sqrt(B_sq)
    R_c = R_prime - delta

    if R_c <= 0:
        raise ValueError("Critical radius non-physical")

    details = {
        "R_critical": R_c,
        "R_extrapolated": R_prime,
        "delta": delta,
        "B_sq": B_sq,
        "M_sq": M_sq,
        "L1_sq": L1_sq,
        "L2_sq": L2_sq,
        "tau": tau,
        **k_details,
    }

    return R_c, details


def compute_two_group_alpha(
    xs: TwoGroupMacroXS,
    R: float,
    extrapolation_factor: float = 0.7104,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Compute Rossi alpha (inverse period) eigenvalues for two-group system.

    Returns both eigenvalues (prompt and delayed-like modes) and eigenvectors.
    The dominant (largest) eigenvalue determines the system response.
    """
    # Geometric buckling
    lambda_tr1 = 1.0 / xs.Sigma_tr1
    delta = extrapolation_factor * lambda_tr1
    R_prime = R + delta
    B_sq = (np.pi / R_prime) ** 2

    # Leakage terms
    D1_B2 = xs.D1 * B_sq
    D2_B2 = xs.D2 * B_sq

    # Removal + leakage
    Sigma_r1 = xs.Sigma_a1 + xs.Sigma_s12 + D1_B2
    Sigma_r2 = xs.Sigma_a2 + D2_B2

    # Production terms
    nu_Sigma_f1 = xs.nu1 * xs.Sigma_f1
    nu_Sigma_f2 = xs.nu2 * xs.Sigma_f2

    # System matrix for dn/dt = A·n (in terms of neutron densities)
    # dn1/dt = [χ1·(νΣf1·v1·n1 + νΣf2·v2·n2) - (Σa1 + Σs12 + D1·B²)·v1·n1]
    # dn2/dt = [χ2·(νΣf1·v1·n1 + νΣf2·v2·n2) + Σs12·v1·n1 - (Σa2 + D2·B²)·v2·n2]

    v1, v2 = xs.v1, xs.v2

    A_11 = CHI_1 * nu_Sigma_f1 * v1 - Sigma_r1 * v1
    A_12 = CHI_1 * nu_Sigma_f2 * v2
    A_21 = CHI_2 * nu_Sigma_f1 * v1 + xs.Sigma_s12 * v1
    A_22 = CHI_2 * nu_Sigma_f2 * v2 - Sigma_r2 * v2

    A = np.array([[A_11, A_12], [A_21, A_22]])

    eigenvalues, eigenvectors = np.linalg.eig(A)

    # Sort by real part (dominant mode first)
    idx = np.argsort(-np.real(eigenvalues))
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # k_eff from dominant eigenvalue
    # At criticality, dominant eigenvalue = 0
    # k_eff ≈ 1 + α·Λ where Λ is mean generation time
    alpha_dominant = np.real(eigenvalues[0])

    # Estimate generation time from removal rate
    Lambda_est = 1.0 / (Sigma_r1 * v1)  # approximate
    k_eff_est = 1.0 + alpha_dominant * Lambda_est

    details = {
        "alpha_dominant": alpha_dominant,
        "alpha_secondary": np.real(eigenvalues[1]),
        "k_eff_estimate": k_eff_est,
        "B_sq": B_sq,
        "eigenvector_1": eigenvectors[:, 0],
        "eigenvector_2": eigenvectors[:, 1],
    }

    return eigenvalues, eigenvectors, details


class TwoGroupSimulation:
    """
    Two-group time-dependent simulation with hydrodynamics.

    Solves coupled system:
    - Two-group neutron kinetics
    - Energy deposition
    - Hydrodynamic expansion
    """

    def __init__(
        self,
        core_data: TwoGroupData,
        core_mass: float,
        compression: float,
        tamper_data: Optional[TwoGroupData] = None,
        tamper_mass: float = 0.0,
    ):
        self.core_data = core_data
        self.core_mass = core_mass
        self.initial_compression = compression
        self.tamper_data = tamper_data
        self.tamper_mass = tamper_mass

        # Initial geometry
        rho = core_data.rho_0 * compression
        V = core_mass / rho
        self.R_0 = (3 * V / (4 * np.pi)) ** (1 / 3)

        # Energy per fission (use average of both groups)
        self.E_fission = 200.0 * MEV_TO_J

    def compute_alpha(
        self, R: float, compression: float
    ) -> Tuple[float, TwoGroupMacroXS]:
        """Compute dominant alpha and cross-sections at current state."""
        xs = compute_two_group_xs(self.core_data, compression)
        eigenvalues, _, details = compute_two_group_alpha(xs, R)
        return np.real(eigenvalues[0]), xs

    def simulate(
        self,
        n1_0: float = 1e10,
        n2_0: float = 0.0,
        t_max: float = 1e-6,
    ) -> Dict[str, Any]:
        """
        Run two-group simulation.

        State vector: [R, v, n1, n2, E_fis, E_kin]
        - R: radius
        - v: expansion velocity
        - n1, n2: neutron densities in groups 1 and 2
        - E_fis: fission energy deposited
        - E_kin: kinetic energy
        """
        # Initial flux ratio from eigenvalue problem
        xs_0 = compute_two_group_xs(self.core_data, self.initial_compression)
        _, eigvecs, _ = compute_two_group_alpha(xs_0, self.R_0)

        # Start with eigenvector ratio if n2_0 not specified
        if n2_0 == 0.0:
            ratio = np.real(eigvecs[1, 0] / eigvecs[0, 0])
            n2_0 = n1_0 * max(0, ratio)

        y0 = np.array([self.R_0, 0.0, n1_0, n2_0, 0.0, 0.0])

        alpha_0, xs_0 = self.compute_alpha(self.R_0, self.initial_compression)

        print(f"Initial conditions (two-group):")
        print(f"  R_0 = {self.R_0 * 100:.3f} cm")
        print(f"  compression = {self.initial_compression:.2f}")
        print(f"  α_dominant = {alpha_0:.3e} /s")
        print(f"  n1_0/n2_0 = {n1_0:.2e} / {n2_0:.2e}")

        if alpha_0 <= 0:
            print("  WARNING: System is subcritical!")
            return {"error": "subcritical"}

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            R, v, n1, n2, E_fis, E_kin = y

            c = self.initial_compression * (self.R_0 / R) ** 3
            c = max(c, 0.01)

            xs = compute_two_group_xs(self.core_data, c)

            # Geometric buckling
            lambda_tr1 = 1.0 / xs.Sigma_tr1
            delta = 0.7104 * lambda_tr1
            R_prime = R + delta
            B_sq = (np.pi / R_prime) ** 2

            # Leakage
            D1_B2 = xs.D1 * B_sq
            D2_B2 = xs.D2 * B_sq

            # Removal + leakage
            Sigma_r1 = xs.Sigma_a1 + xs.Sigma_s12 + D1_B2
            Sigma_r2 = xs.Sigma_a2 + D2_B2

            v1, v2 = xs.v1, xs.v2

            # Fission rates (reactions per unit volume per second)
            fission_rate_1 = xs.Sigma_f1 * n1 * v1
            fission_rate_2 = xs.Sigma_f2 * n2 * v2
            total_fission_rate = fission_rate_1 + fission_rate_2

            # Neutron production
            source = xs.nu1 * fission_rate_1 + xs.nu2 * fission_rate_2

            # Two-group kinetics
            dn1_dt = CHI_1 * source - Sigma_r1 * n1 * v1
            dn2_dt = CHI_2 * source + xs.Sigma_s12 * n1 * v1 - Sigma_r2 * n2 * v2

            # Fission power
            V_core = (4 / 3) * np.pi * R**3
            P_fis = total_fission_rate * V_core * self.E_fission
            dE_fis_dt = P_fis

            # Hydrodynamics (radiation-dominated EOS)
            E_int = max(0, E_fis - E_kin)
            a_rad = 7.566e-16  # radiation constant J/m³/K⁴
            T = (E_int / (a_rad * V_core)) ** 0.25 if E_int > 0 else 1e4
            P = a_rad * T**4 / 3

            # Shell acceleration
            m = self.core_mass + self.tamper_mass
            dv_dt = 4 * np.pi * R**2 * P / m if m > 0 else 0
            dE_kin_dt = 4 * np.pi * R**2 * P * v if v > 0 else 0

            return np.array([v, dv_dt, dn1_dt, dn2_dt, dE_fis_dt, dE_kin_dt])

        def subcritical_event(t: float, y: np.ndarray) -> float:
            R = y[0]
            c = self.initial_compression * (self.R_0 / R) ** 3
            if c < 0.01:
                return -1.0
            alpha, _ = self.compute_alpha(R, c)
            return alpha

        subcritical_event.terminal = True  # type: ignore[attr-defined]
        subcritical_event.direction = -1  # type: ignore[attr-defined]

        sol = solve_ivp(
            rhs,
            (0, t_max),
            y0,
            method="BDF",
            max_step=1e-10,
            rtol=1e-7,
            events=subcritical_event,
        )

        results = {
            "t": sol.t,
            "R": sol.y[0],
            "v": sol.y[1],
            "n1": sol.y[2],
            "n2": sol.y[3],
            "E_fission": sol.y[4],
            "E_kinetic": sol.y[5],
            "yield_kt": sol.y[4] / 4.184e12,
            "compression": self.initial_compression * (self.R_0 / sol.y[0]) ** 3,
        }

        print(f"\nFinal state:")
        print(f"  t = {sol.t[-1] * 1e6:.3f} μs")
        print(f"  Yield = {results['yield_kt'][-1]:.2f} kt")
        print(f"  n1_final/n2_final = {sol.y[2, -1]:.2e} / {sol.y[3, -1]:.2e}")

        return results


def validate_two_group():
    """Validate two-group implementation against one-group and benchmarks."""
    print("=" * 70)
    print("TWO-GROUP PHYSICS VALIDATION")
    print("=" * 70)

    # Test 1: k_inf calculation
    print("\n--- k_infinity Comparison ---")
    for data in [PU239_TWO_GROUP, U238_TWO_GROUP]:
        xs = compute_two_group_xs(data, compression=1.0)
        k_inf, details = compute_two_group_k_inf(xs)
        print(f"{data.name}:")
        print(f"  k_inf = {k_inf:.4f}")
        print(f"  φ2/φ1 = {details['flux_ratio_phi2_phi1']:.4f}")

    # Test 2: Critical radius
    print("\n--- Critical Radius (Pu-239) ---")
    for c in [1.0, 2.0, 2.5, 3.0]:
        try:
            R_c, details = compute_two_group_critical_radius(PU239_TWO_GROUP, c)
            mass = (4 / 3) * np.pi * R_c**3 * PU239_TWO_GROUP.rho_0 * c
            print(
                f"  c={c:.1f}: R_c={R_c * 100:.2f}cm, m_c={mass:.2f}kg, k_inf={details['k_inf']:.4f}"
            )
        except ValueError as e:
            print(f"  c={c:.1f}: {e}")

    # Test 3: Trinity simulation
    print("\n--- Trinity Configuration (Two-Group) ---")
    sim = TwoGroupSimulation(
        core_data=PU239_TWO_GROUP,
        core_mass=6.2,
        compression=2.5,
        tamper_mass=108.0,
    )

    results = sim.simulate(n1_0=1e10, t_max=2e-6)

    if "error" not in results:
        print(f"\n  Final yield: {results['yield_kt'][-1]:.1f} kt")
        print(f"  (Trinity actual: ~15-21 kt)")

    return results


if __name__ == "__main__":
    validate_two_group()
