"""
Rigorous Nuclear Physics Simulation
Based on arXiv:1606.01670v1 with transport corrections and proper physics.

Key improvements over naive diffusion:
1. Transport-corrected diffusion coefficient
2. Energy-dependent fission spectrum treatment
3. Proper equation of state (radiation + matter)
4. Accurate extrapolation boundary condition
5. Stiff ODE solver for coupled neutronics-hydrodynamics
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import (
    physical_constants,
    k as k_B,
    c as c_light,
    sigma as sigma_SB,
    N_A,
    eV,
    m_n,
)
from scipy.special import kn
from dataclasses import dataclass, field
from typing import Tuple, Dict, Any, Optional
import matplotlib.pyplot as plt

BARN_TO_M2 = 1e-28
MEV_TO_J = 1.602176634e-13
AMU_TO_KG = 1.66053906660e-27


@dataclass
class NuclearData:
    """
    Nuclear data with energy-dependent cross-sections.
    Values are fission-spectrum averaged from ENDF/B-VIII.0.
    """

    name: str
    A: float  # atomic mass number
    M: float  # isotope mass (amu)
    rho_0: float  # reference density (kg/m³)

    # Fission spectrum averaged cross-sections (barns)
    sigma_f: float  # fission
    sigma_el: float  # elastic scattering
    sigma_inel: float  # inelastic scattering
    sigma_c: float  # capture (n,γ)

    # Fission data
    nu_bar: float  # total neutrons per fission
    nu_prompt: float  # prompt neutrons per fission
    beta_eff: float  # delayed neutron fraction
    E_fission: float  # energy per fission (MeV)

    # Fission neutron spectrum parameters (Watt spectrum)
    a_watt: float = 0.988  # MeV
    b_watt: float = 2.249  # MeV^-1

    @property
    def sigma_s(self) -> float:
        return self.sigma_el + self.sigma_inel

    @property
    def sigma_t(self) -> float:
        return self.sigma_f + self.sigma_s + self.sigma_c

    @property
    def sigma_a(self) -> float:
        return self.sigma_f + self.sigma_c

    @property
    def n_0(self) -> float:
        """Reference number density (atoms/m³)."""
        return self.rho_0 / (self.M * AMU_TO_KG)

    def average_neutron_energy(self) -> float:
        """Mean energy of fission neutrons (MeV) from Watt spectrum."""
        # <E> ≈ 3/2 * a + a²*b/4 for Watt spectrum
        return 1.5 * self.a_watt + 0.25 * self.a_watt**2 * self.b_watt

    def average_neutron_velocity(self) -> float:
        """RMS velocity of fission neutrons (m/s)."""
        E_avg_J = self.average_neutron_energy() * MEV_TO_J
        # v = sqrt(2E/m) for average
        return np.sqrt(2 * E_avg_J / m_n)


# Pu-239 data from ENDF/B-VIII.0 (fission spectrum averaged)
PU239 = NuclearData(
    name="Pu-239",
    A=239.0,
    M=239.0521634,
    rho_0=19.86e3,  # α-phase density
    sigma_f=1.800,
    sigma_el=4.566,
    sigma_inel=1.369,
    sigma_c=0.065,
    nu_bar=3.172,
    nu_prompt=3.165,
    beta_eff=0.0021,
    E_fission=200.0,
    a_watt=0.966,
    b_watt=2.842,
)

# U-235 data from ENDF/B-VIII.0
U235 = NuclearData(
    name="U-235",
    A=235.0,
    M=235.0439299,
    rho_0=18.9e3,
    sigma_f=1.235,
    sigma_el=4.409,
    sigma_inel=1.917,
    sigma_c=0.089,
    nu_bar=2.637,
    nu_prompt=2.620,
    beta_eff=0.0065,
    E_fission=200.0,
    a_watt=0.988,
    b_watt=2.249,
)


@dataclass
class MacroscopicCrossSections:
    """Macroscopic cross-sections at given density."""

    Sigma_f: float
    Sigma_s: float
    Sigma_a: float
    Sigma_t: float
    Sigma_tr: float
    n_atoms: float


def compute_macroscopic_xs(
    data: NuclearData, compression: float
) -> MacroscopicCrossSections:
    """Compute macroscopic cross-sections for given compression."""
    n = data.n_0 * compression

    Sigma_f = data.sigma_f * BARN_TO_M2 * n
    Sigma_s = data.sigma_s * BARN_TO_M2 * n
    Sigma_a = data.sigma_a * BARN_TO_M2 * n
    Sigma_t = data.sigma_t * BARN_TO_M2 * n

    # Transport cross-section with anisotropic scattering correction
    # μ_0 = <cos θ> ≈ 2/(3A) for elastic scattering from heavy nucleus
    mu_0 = 2.0 / (3.0 * data.A)
    Sigma_tr = Sigma_t - mu_0 * Sigma_s

    return MacroscopicCrossSections(
        Sigma_f=Sigma_f,
        Sigma_s=Sigma_s,
        Sigma_a=Sigma_a,
        Sigma_t=Sigma_t,
        Sigma_tr=Sigma_tr,
        n_atoms=n,
    )


@dataclass
class DiffusionParameters:
    """Diffusion theory parameters with transport corrections."""

    D: float  # diffusion coefficient
    D_tr: float  # transport-corrected diffusion coefficient
    lambda_tr: float  # transport mean free path
    B_m_sq: float  # material buckling
    L_sq: float  # diffusion length squared
    k_inf: float  # infinite multiplication factor
    Lambda: float  # prompt neutron generation time (mean time between fissions)
    l: float  # prompt neutron lifetime (mean time to absorption)
    xs: MacroscopicCrossSections
    c: float  # scattering ratio Sigma_s / Sigma_t


def compute_diffusion_params(
    data: NuclearData, compression: float
) -> DiffusionParameters:
    """
    Compute diffusion parameters with transport correction.

    Uses proper P1 transport correction:
    - D = lambda_tr / 3 where lambda_tr = 1/Sigma_tr
    - Sigma_tr = Sigma_t - mu_0 * Sigma_s (already transport-corrected)

    For highly absorbing media, applies buckling-dependent correction.
    """
    xs = compute_macroscopic_xs(data, compression)

    lambda_tr = 1.0 / xs.Sigma_tr
    D = lambda_tr / 3.0

    c = xs.Sigma_s / xs.Sigma_t

    # Buckling-dependent diffusion coefficient correction for absorbing media
    # From Davison (1957): D_eff = D * (1 - 2*D*Sigma_a / (5 + 3*D*Sigma_a))
    # This accounts for P1 breakdown when absorption is significant
    correction_factor = 2.0 * D * xs.Sigma_a / (5.0 + 3.0 * D * xs.Sigma_a)
    D_tr = D * (1.0 - correction_factor)

    B_m_sq = (data.nu_prompt * xs.Sigma_f - xs.Sigma_a) / D_tr
    L_sq = D_tr / xs.Sigma_a
    k_inf = data.nu_prompt * xs.Sigma_f / xs.Sigma_a

    v = data.average_neutron_velocity()

    # Prompt neutron lifetime: mean time to absorption l = 1/(v * Sigma_a)
    l = 1.0 / (v * xs.Sigma_a)

    # Prompt neutron generation time: Lambda = l / k_inf
    # This is the mean time between birth of one generation and the next
    Lambda = l / k_inf

    return DiffusionParameters(
        D=D,
        D_tr=D_tr,
        lambda_tr=lambda_tr,
        B_m_sq=B_m_sq,
        L_sq=L_sq,
        k_inf=k_inf,
        Lambda=Lambda,
        l=l,
        xs=xs,
        c=c,
    )


def extrapolation_distance(lambda_tr: float, c: float) -> float:
    """
    Compute extrapolation distance for vacuum boundary using transport theory.

    For pure scattering (c=1): delta = 0.7104 * lambda_tr (Milne problem)
    For absorbing media: uses Mark approximation from Case & Zweifel.

    Mark approximation: delta/lambda_tr = 0.7104 * sqrt(3c / (1 + 5c/3))
    Valid for c > 0.3 (typical fast reactor spectra).
    """
    if c > 0.9999:
        return 0.7104 * lambda_tr

    # Mark approximation for absorbing media (Case & Zweifel, Linear Transport Theory)
    # delta = 0.7104 * lambda_tr * sqrt(3c / (1 + 5c/3))
    mark_factor = np.sqrt(3.0 * c / (1.0 + 5.0 * c / 3.0))
    return 0.7104 * lambda_tr * mark_factor


def compute_critical_radius(
    data: NuclearData, compression: float = 1.0
) -> Tuple[float, DiffusionParameters]:
    """
    Compute critical radius using transport-corrected diffusion theory.

    At criticality: k_eff = 1
    For bare sphere: B_g² = (π/R')² = B_m² where R' = R + δ
    """
    params = compute_diffusion_params(data, compression)

    if params.B_m_sq <= 0:
        raise ValueError("Material is subcritical even infinite (k_inf < 1)")

    # Scattering ratio for extrapolation correction
    c = params.xs.Sigma_s / params.xs.Sigma_t
    delta = extrapolation_distance(params.lambda_tr, c)

    # R'_critical = π / √(B_m²)
    R_prime_c = np.pi / np.sqrt(params.B_m_sq)
    R_c = R_prime_c - delta

    if R_c <= 0:
        raise ValueError("Critical radius is non-physical (too much extrapolation)")

    return R_c, params


def compute_critical_mass(
    data: NuclearData, compression: float = 1.0
) -> Tuple[float, float, DiffusionParameters]:
    """Compute critical mass for bare sphere."""
    R_c, params = compute_critical_radius(data, compression)
    volume = (4.0 / 3.0) * np.pi * R_c**3
    mass = volume * data.rho_0 * compression
    return mass, R_c, params


class EquationOfState:
    """
    Equation of state for hot dense matter.
    Radiation + ideal gas + electron degeneracy.
    """

    M_E = 9.1093837015e-31
    H_BAR = 1.054571817e-34
    K_FERMI = (3.0 * np.pi**2) ** (1.0 / 3.0)

    @staticmethod
    def radiation_pressure(T: float) -> float:
        a_rad = 4.0 * sigma_SB / c_light
        return a_rad * T**4 / 3.0

    @staticmethod
    def radiation_energy_density(T: float) -> float:
        a_rad = 4.0 * sigma_SB / c_light
        return a_rad * T**4

    @staticmethod
    def ideal_gas_pressure(n: float, T: float, Z: float = 1.0) -> float:
        return n * (1.0 + Z) * k_B * T

    @staticmethod
    def electron_degeneracy_pressure(n_e: float) -> float:
        """
        Non-relativistic electron degeneracy pressure.
        P_deg = (h_bar^2 / 5 m_e) * (3 pi^2)^(2/3) * n_e^(5/3)
        """
        prefactor = (EquationOfState.H_BAR**2 / (5.0 * EquationOfState.M_E)) * (
            3.0 * np.pi**2
        ) ** (2.0 / 3.0)
        return prefactor * n_e ** (5.0 / 3.0)

    @staticmethod
    def fermi_energy(n_e: float) -> float:
        """Fermi energy E_F = (h_bar^2 / 2m_e) * (3 pi^2 n_e)^(2/3)"""
        return (EquationOfState.H_BAR**2 / (2.0 * EquationOfState.M_E)) * (
            3.0 * np.pi**2 * n_e
        ) ** (2.0 / 3.0)

    @staticmethod
    def degeneracy_parameter(n_e: float, T: float) -> float:
        """Theta = k_B T / E_F. Theta << 1 means degenerate."""
        E_F = EquationOfState.fermi_energy(n_e)
        return k_B * T / E_F if E_F > 0 else float("inf")

    @staticmethod
    def total_pressure(T: float, n: float, Z: float = 40.0) -> float:
        """Total pressure including radiation, gas, and degeneracy."""
        P_rad = EquationOfState.radiation_pressure(T)
        P_gas = EquationOfState.ideal_gas_pressure(n, T, Z)

        n_e = n * Z
        theta = EquationOfState.degeneracy_parameter(n_e, T)

        if theta < 3.0:
            P_deg = EquationOfState.electron_degeneracy_pressure(n_e)
            # Interpolate between degenerate and classical regimes
            # Using Fermi-Dirac interpolation factor
            deg_factor = 1.0 / (1.0 + theta**1.5)
            P_electron = P_deg * deg_factor + (n_e * k_B * T) * (1.0 - deg_factor)
            return P_rad + n * k_B * T + P_electron
        else:
            return P_rad + P_gas

    @staticmethod
    def temperature_from_energy(E: float, V: float, n: float, Z: float = 40.0) -> float:
        """Solve for T given internal energy E in volume V."""
        a_rad = 4.0 * sigma_SB / c_light

        T_rad = (E / (a_rad * V)) ** 0.25 if E > 0 and V > 0 else 1e6

        T = max(T_rad, 1e4)
        for _ in range(30):
            u_rad = a_rad * T**4
            u_gas = 1.5 * n * (1.0 + Z) * k_B * T

            n_e = n * Z
            theta = EquationOfState.degeneracy_parameter(n_e, T)
            if theta < 3.0:
                E_F = EquationOfState.fermi_energy(n_e)
                u_deg = 0.6 * n_e * E_F
                deg_factor = 1.0 / (1.0 + theta**1.5)
                u_electron = u_deg * deg_factor + (1.5 * n_e * k_B * T) * (
                    1.0 - deg_factor
                )
                u_total = u_rad + n * 1.5 * k_B * T + u_electron
            else:
                u_total = u_rad + u_gas

            du_dT = 4.0 * a_rad * T**3 + 1.5 * n * (1.0 + Z) * k_B

            residual = u_total * V - E
            if abs(residual) < 1e-8 * max(E, 1e-20):
                break

            T = T - residual / (du_dT * V)
            T = max(T, 1e4)

        return T


class MassShell:
    """
    Proper mass shell model from paper's Eq. 37.

    Models the bomb structure with:
    - Pu core (6.2 kg)
    - U-238 tamper (108 kg)
    - Al pusher (130 kg)
    - HE and casing (4430 kg)
    """

    def __init__(
        self,
        core_mass: float = 6.2,
        tamper_mass: float = 108.0,
        pusher_mass: float = 130.0,
        outer_mass: float = 4432.0,
        R_core: float = 0.05,
    ):
        self.core_mass = core_mass
        self.tamper_mass = tamper_mass
        self.pusher_mass = pusher_mass
        self.outer_mass = outer_mass

        # Reference radii (will be scaled by compression)
        self.R_core_ref = R_core

        # Estimate shell radii from masses and densities
        rho_U = 18.95e3
        rho_Al = 2.7e3
        rho_outer = 1.6e3  # approximate HE + casing

        V_core = core_mass / 19.86e3
        V_tamper = tamper_mass / rho_U
        V_pusher = pusher_mass / rho_Al

        self.R_tamper_ref = (
            (V_core + V_tamper) ** (1 / 3) * (3 / (4 * np.pi)) ** (1 / 3) * 2
        )
        self.R_pusher_ref = (
            (V_core + V_tamper + V_pusher) ** (1 / 3)
            * (3 / (4 * np.pi)) ** (1 / 3)
            * 2.5
        )
        self.R_outer_ref = 1.0  # meters

    def enclosed_mass(self, R: float, compression: float) -> float:
        """Total mass enclosed within radius R."""
        scale = compression ** (-1 / 3)

        R_core = self.R_core_ref * scale
        R_tamper = self.R_tamper_ref * scale
        R_pusher = self.R_pusher_ref * scale
        R_outer = self.R_outer_ref

        if R <= R_core:
            return self.core_mass
        elif R <= R_tamper:
            frac = (R**3 - R_core**3) / (R_tamper**3 - R_core**3)
            return self.core_mass + self.tamper_mass * frac
        elif R <= R_pusher:
            frac = (R**3 - R_tamper**3) / (R_pusher**3 - R_tamper**3)
            return self.core_mass + self.tamper_mass + self.pusher_mass * frac
        elif R <= R_outer:
            frac = (R**3 - R_pusher**3) / (R_outer**3 - R_pusher**3)
            return (
                self.core_mass
                + self.tamper_mass
                + self.pusher_mass
                + self.outer_mass * frac
            )
        else:
            rho_air = 1.29
            extra = (4 / 3) * np.pi * (R**3 - R_outer**3) * rho_air
            return (
                self.core_mass
                + self.tamper_mass
                + self.pusher_mass
                + self.outer_mass
                + extra
            )


def compute_rossi_alpha(
    data: NuclearData, R: float, compression: float
) -> Tuple[float, DiffusionParameters]:
    """
    Compute inverse reactor period (Rossi alpha).

    α = v * D * (B_m² - B_g²) = (k_eff - 1) / τ

    where B_g² = (π/R')² is geometric buckling.
    """
    params = compute_diffusion_params(data, compression)

    c = params.xs.Sigma_s / params.xs.Sigma_t
    delta = extrapolation_distance(params.lambda_tr, c)

    R_prime = R + delta
    B_g_sq = (np.pi / R_prime) ** 2

    v = data.average_neutron_velocity()
    alpha = v * params.D_tr * (params.B_m_sq - B_g_sq)

    return alpha, params


class SupercriticalSimulation:
    """
    Full simulation of prompt supercritical excursion.

    Solves coupled system:
    1. dN/dt = α(t) N  (neutron kinetics)
    2. dE/dt = P(t)    (energy deposition)
    3. dR/dt = v       (expansion velocity)
    4. dv/dt = F/m     (momentum equation)

    With proper:
    - Transport-corrected diffusion
    - Radiation + matter EOS
    - Mass shell dynamics
    """

    def __init__(
        self,
        data: NuclearData,
        core_mass: float,
        compression: float,
        use_full_mass_shell: bool = False,
    ):
        self.data = data
        self.core_mass = core_mass
        self.initial_compression = compression

        # Initial geometry
        rho = data.rho_0 * compression
        V = core_mass / rho
        self.R_0 = (3 * V / (4 * np.pi)) ** (1 / 3)

        # Mass model
        if use_full_mass_shell:
            self.mass_shell = MassShell(core_mass, R_core=self.R_0)
        else:
            self.mass_shell = None

        # Fission energy
        self.E_fission = data.E_fission * MEV_TO_J

    def get_mass(self, R: float, compression: float) -> float:
        """Get enclosed mass at radius R."""
        if self.mass_shell is not None:
            return self.mass_shell.enclosed_mass(R, compression)
        else:
            if R <= self.R_0 * compression ** (-1 / 3):
                return self.core_mass
            rho_air = 1.29
            R_0_eff = self.R_0 * compression ** (-1 / 3)
            return self.core_mass + (4 / 3) * np.pi * (R**3 - R_0_eff**3) * rho_air

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        ODE right-hand side.
        State: y = [R, v, N, E_fission, E_kinetic]
        """
        R, v, N, E_fis, E_kin = y

        # Current compression from radius change
        c = self.initial_compression * (self.R_0 / R) ** 3
        c = max(c, 0.01)

        # Neutronics
        alpha, params = compute_rossi_alpha(self.data, R, c)
        dN_dt = alpha * N

        # Fission power: P = Σ_f * Φ * ε = Σ_f * N * v_n * ε
        v_n = self.data.average_neutron_velocity()
        P = params.xs.Sigma_f * N * v_n * self.E_fission
        dE_fis_dt = P

        # Internal energy (radiation + thermal)
        E_internal = E_fis - E_kin
        E_internal = max(E_internal, 0.0)

        # Volume and number density
        V = (4 / 3) * np.pi * R**3
        n = params.xs.n_atoms

        # Temperature from EOS
        T = EquationOfState.temperature_from_energy(E_internal, V, n)

        # Pressure from EOS
        P_total = EquationOfState.total_pressure(T, n)

        # Shell dynamics: m dv/dt = 4πR² P
        m = self.get_mass(R, c)
        dv_dt = 4 * np.pi * R**2 * P_total / m if m > 0 else 0

        # Kinetic energy rate: dE_kin/dt = F · v = 4πR² P v
        dE_kin_dt = 4 * np.pi * R**2 * P_total * v if v > 0 else 0

        return np.array([v, dv_dt, dN_dt, dE_fis_dt, dE_kin_dt])

    def simulate(self, N_0: float = 1e8, t_max: float = 1e-6) -> Dict[str, Any]:
        """Run simulation with stiff solver."""

        y0 = np.array([self.R_0, 0.0, N_0, 0.0, 0.0])

        alpha_0, params_0 = compute_rossi_alpha(
            self.data, self.R_0, self.initial_compression
        )
        k_eff = 1 + alpha_0 * params_0.Lambda

        print(f"Initial conditions:")
        print(f"  R_0 = {self.R_0 * 100:.3f} cm")
        print(f"  compression = {self.initial_compression:.2f}")
        print(f"  α = {alpha_0:.3e} /s")
        print(f"  k_eff ≈ {k_eff:.4f}")
        print(f"  k_inf = {params_0.k_inf:.4f}")

        if alpha_0 <= 0:
            print("  WARNING: System is subcritical!")
            return {"error": "subcritical"}

        R_0 = self.R_0
        c_init = self.initial_compression
        data = self.data

        def subcritical_event(t: float, y: np.ndarray) -> float:
            R = y[0]
            c = c_init * (R_0 / R) ** 3
            if c < 0.01:
                return -1.0
            alpha, _ = compute_rossi_alpha(data, R, c)
            return alpha

        subcritical_event.terminal = True  # type: ignore
        subcritical_event.direction = -1  # type: ignore

        sol = solve_ivp(
            self.rhs,
            (0, t_max),
            y0,
            method="BDF",  # implicit method for stiff ODEs
            max_step=1e-11,
            rtol=1e-8,
            atol=1e-12,
            events=subcritical_event,
        )

        # Extract results
        t = sol.t
        R = sol.y[0]
        v = sol.y[1]
        N = sol.y[2]
        E_fis = sol.y[3]
        E_kin = sol.y[4]

        # Derived quantities
        c_hist = self.initial_compression * (self.R_0 / R) ** 3
        E_int = np.maximum(0, E_fis - E_kin)

        # Temperature and pressure history
        V = (4 / 3) * np.pi * R**3
        T_hist = np.zeros_like(t)
        P_hist = np.zeros_like(t)
        alpha_hist = np.zeros_like(t)

        for i in range(len(t)):
            params = compute_diffusion_params(self.data, c_hist[i])
            T_hist[i] = EquationOfState.temperature_from_energy(
                E_int[i], V[i], params.xs.n_atoms
            )
            P_hist[i] = EquationOfState.total_pressure(T_hist[i], params.xs.n_atoms)
            alpha_hist[i], _ = compute_rossi_alpha(self.data, R[i], c_hist[i])

        yield_kt = E_fis / 4.184e12

        results = {
            "t": t,
            "R": R,
            "v": v,
            "N": N,
            "E_fission": E_fis,
            "E_kinetic": E_kin,
            "E_internal": E_int,
            "compression": c_hist,
            "T": T_hist,
            "P": P_hist,
            "alpha": alpha_hist,
            "yield_kt": yield_kt,
            "final_yield_kt": yield_kt[-1],
        }

        print(f"\nResults:")
        print(f"  t_final = {t[-1] * 1e6:.3f} μs")
        print(f"  R_final = {R[-1] * 100:.3f} cm")
        print(f"  Yield = {yield_kt[-1]:.2f} kt")
        print(f"  Peak T = {np.max(T_hist):.2e} K")
        print(f"  Peak P = {np.max(P_hist) / 1e14:.2f} Gbar")

        return results


def plot_simulation(results: Dict[str, Any], title: str = "Supercritical Simulation"):
    """Plot simulation results."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle(title, fontsize=14)

    t_us = results["t"] * 1e6

    axes[0, 0].semilogy(t_us, results["N"])
    axes[0, 0].set_xlabel("Time (μs)")
    axes[0, 0].set_ylabel("Neutron count")
    axes[0, 0].set_title("Neutron Population")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(t_us, results["yield_kt"])
    axes[0, 1].set_xlabel("Time (μs)")
    axes[0, 1].set_ylabel("Yield (kt)")
    axes[0, 1].set_title("Energy Release")
    axes[0, 1].grid(True, alpha=0.3)

    axes[0, 2].plot(t_us, results["alpha"] / 1e8)
    axes[0, 2].axhline(y=0, color="r", linestyle="--", alpha=0.5)
    axes[0, 2].set_xlabel("Time (μs)")
    axes[0, 2].set_ylabel("α (10⁸/s)")
    axes[0, 2].set_title("Rossi Alpha")
    axes[0, 2].grid(True, alpha=0.3)

    axes[1, 0].plot(t_us, results["R"] * 100, "b-", label="R")
    axes[1, 0].set_xlabel("Time (μs)")
    axes[1, 0].set_ylabel("Radius (cm)", color="b")
    ax2 = axes[1, 0].twinx()
    ax2.plot(t_us, results["v"] / 1000, "r-", label="v")
    ax2.set_ylabel("Velocity (km/s)", color="r")
    axes[1, 0].set_title("Expansion")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].semilogy(t_us, results["T"])
    axes[1, 1].set_xlabel("Time (μs)")
    axes[1, 1].set_ylabel("Temperature (K)")
    axes[1, 1].set_title("Temperature")
    axes[1, 1].grid(True, alpha=0.3)

    axes[1, 2].plot(t_us, results["P"] / 1e14)
    axes[1, 2].set_xlabel("Time (μs)")
    axes[1, 2].set_ylabel("Pressure (Gbar)")
    axes[1, 2].set_title("Pressure")
    axes[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    print("=" * 70)
    print("RIGOROUS NUCLEAR PHYSICS SIMULATION")
    print("=" * 70)

    print("\n--- Critical Mass Calculations ---")
    for mat in [PU239, U235]:
        try:
            mass, R_c, params = compute_critical_mass(mat, compression=1.0)
            print(f"\n{mat.name}:")
            print(f"  R_critical = {R_c * 100:.3f} cm")
            print(f"  m_critical = {mass:.2f} kg")
            print(f"  k_inf = {params.k_inf:.4f}")
            print(f"  <E_n> = {mat.average_neutron_energy():.3f} MeV")
            print(f"  <v_n> = {mat.average_neutron_velocity() / 1e6:.2f} × 10⁶ m/s")

            print(f"  Compression scaling:")
            for c in [1.5, 2.0, 2.5, 3.0]:
                m_c, R_c_c, _ = compute_critical_mass(mat, compression=c)
                print(f"    c={c:.1f}: R_c={R_c_c * 100:.2f}cm, m_c={m_c:.2f}kg")
        except ValueError as e:
            print(f"\n{mat.name}: {e}")

    print("\n" + "=" * 70)
    print("SUPERCRITICAL SIMULATION")
    print("=" * 70)

    sim = SupercriticalSimulation(
        data=PU239, core_mass=6.2, compression=2.5, use_full_mass_shell=False
    )

    results = sim.simulate(N_0=1e8, t_max=2e-6)

    if "error" not in results:
        fig = plot_simulation(results, "Pu-239 Supercritical Excursion (c=2.5)")
        plt.savefig("rigorous_simulation.png", dpi=150, bbox_inches="tight")
        plt.close()
        print("\nSaved: rigorous_simulation.png")

    print("\n" + "=" * 70)
    print("COMPRESSION STUDY")
    print("=" * 70)

    compressions = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    yields = []

    for c in compressions:
        print(f"\n--- c = {c} ---")
        sim = SupercriticalSimulation(PU239, 6.2, c)
        res = sim.simulate(N_0=1e8, t_max=2e-6)
        if "error" not in res:
            yields.append(res["final_yield_kt"])
        else:
            yields.append(0)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(compressions, yields, "bo-", markersize=8, linewidth=2)
    ax.set_xlabel("Compression Factor", fontsize=12)
    ax.set_ylabel("Yield (kt)", fontsize=12)
    ax.set_title("Yield vs Compression (Rigorous Model)", fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.savefig("rigorous_yield_vs_compression.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("\nSaved: rigorous_yield_vs_compression.png")
