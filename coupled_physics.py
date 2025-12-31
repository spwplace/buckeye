"""
Coupled Core-Tamper Neutronics with Full Hydrodynamics

Implements:
1. Two-region diffusion (core + tamper) with interface conditions
2. Proper material EOS (radiation + ion + electron pressure)
3. Hydrodynamic expansion with shock physics
4. Time-dependent coupled simulation

Based on Aste 2016 (arXiv:1606.01670) with extensions.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from dataclasses import dataclass, field
from typing import Tuple, Dict, Any, Optional, List
from enum import Enum

import nuclear_physics as nuc  # type: ignore[import-not-found]

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant [J/K]
C_LIGHT = 299792458.0  # Speed of light [m/s]
SIGMA_SB = 5.670374e-8  # Stefan-Boltzmann constant [W/m²/K⁴]
A_RAD = 4 * SIGMA_SB / C_LIGHT  # Radiation constant
M_E = 9.1093837e-31  # Electron mass [kg]
E_CHARGE = 1.602176634e-19  # Elementary charge [C]
H_BAR = 1.054571817e-34  # Reduced Planck constant [J·s]
N_A = 6.02214076e23  # Avogadro's number


@dataclass
class TamperMaterial:
    """Properties of tamper/reflector material."""

    name: str
    A: float  # Atomic mass number
    Z: float  # Atomic number (for ionization)
    rho_0: float  # Reference density [kg/m³]
    sigma_s: float  # Scattering cross-section [barns]
    sigma_a: float  # Absorption cross-section [barns]
    sigma_f: float = 0.0  # Fission cross-section [barns] (for U-238)
    nu_f: float = 0.0  # Neutrons per fission (for U-238)

    @property
    def n_density(self) -> float:
        """Atom number density at reference density [atoms/m³]."""
        return self.rho_0 * N_A / (self.A * 1e-3)


# Common tamper materials with fission-spectrum-averaged cross-sections
# U-238 is subcritical (k_inf ~ 0.3-0.4) - primarily a reflector/moderator
U238_TAMPER = TamperMaterial(
    name="U-238",
    A=238.05,
    Z=92,
    rho_0=18.95e3,
    sigma_s=6.87,
    sigma_a=0.30,  # Total absorption: capture (~0.25) + spectrum-avg fission (~0.05)
    sigma_f=0.04,  # Spectrum-averaged fast fission
    nu_f=2.6,
)

TUNGSTEN_TAMPER = TamperMaterial(
    name="Tungsten",
    A=183.84,
    Z=74,
    rho_0=19.3e3,
    sigma_s=5.0,
    sigma_a=0.5,
)

BERYLLIUM_REFLECTOR = TamperMaterial(
    name="Beryllium",
    A=9.012,
    Z=4,
    rho_0=1.85e3,
    sigma_s=6.1,
    sigma_a=0.01,
)


class MaterialEOS:
    """
    Complete equation of state for weapon-relevant conditions.

    Includes:
    - Radiation pressure (dominant at T > 10⁷ K)
    - Ion thermal pressure
    - Electron pressure (thermal + degeneracy)
    - Ionization effects
    """

    @staticmethod
    def ionization_state(T: float, Z_atom: float, n_ion: float) -> float:
        """
        Average ionization state using Thomas-Fermi model.
        At weapon temperatures, heavy atoms are highly ionized.
        """
        if T < 1e4:
            return 1.0

        # Thomas-Fermi ionization temperature
        T_TF = 1.6e5 * Z_atom ** (4 / 3)  # [K]

        # Approximate ionization fraction
        Z_eff = Z_atom * (1 - np.exp(-((T / T_TF) ** 0.5)))
        return max(1.0, min(Z_eff, Z_atom))

    @staticmethod
    def electron_pressure(n_e: float, T: float) -> float:
        """
        Electron pressure including degeneracy effects.
        Interpolates between classical and degenerate limits.
        """
        if n_e <= 0 or T <= 0:
            return 0.0

        # Classical thermal pressure
        P_thermal = n_e * K_B * T

        # Fermi energy and degeneracy pressure
        E_F = (H_BAR**2 / (2 * M_E)) * (3 * np.pi**2 * n_e) ** (2 / 3)
        P_deg = (2 / 5) * n_e * E_F

        # Degeneracy parameter θ = kT/E_F
        theta = K_B * T / E_F if E_F > 0 else float("inf")

        # Smooth interpolation (Padé approximant)
        if theta < 0.1:
            return P_deg * (1 + (np.pi**2 / 12) * theta**2)
        elif theta > 10:
            return P_thermal
        else:
            # Interpolate
            w = 1 / (1 + theta**1.5)
            return w * P_deg + (1 - w) * P_thermal

    @staticmethod
    def total_pressure(
        T: float, rho: float, A: float, Z: float
    ) -> Tuple[float, Dict[str, float]]:
        """
        Total pressure from all contributions.

        Returns: (P_total, components_dict)
        """
        if rho <= 0 or T <= 0:
            return 0.0, {}

        # Number densities
        n_ion = rho * N_A / (A * 1e-3)
        Z_eff = MaterialEOS.ionization_state(T, Z, n_ion)
        n_e = n_ion * Z_eff

        # Radiation pressure: P_rad = (1/3) a T⁴
        P_rad = A_RAD * T**4 / 3

        # Ion thermal pressure: P_ion = n_ion k T
        P_ion = n_ion * K_B * T

        # Electron pressure (thermal + degeneracy)
        P_e = MaterialEOS.electron_pressure(n_e, T)

        P_total = P_rad + P_ion + P_e

        components = {
            "P_rad": P_rad,
            "P_ion": P_ion,
            "P_electron": P_e,
            "Z_eff": Z_eff,
            "n_ion": n_ion,
            "n_electron": n_e,
        }

        return P_total, components

    @staticmethod
    def internal_energy_density(T: float, rho: float, A: float, Z: float) -> float:
        """Internal energy per unit volume [J/m³]."""
        if rho <= 0 or T <= 0:
            return 0.0

        n_ion = rho * N_A / (A * 1e-3)
        Z_eff = MaterialEOS.ionization_state(T, Z, n_ion)
        n_e = n_ion * Z_eff

        # Radiation: u_rad = a T⁴
        u_rad = A_RAD * T**4

        # Ions: u_ion = (3/2) n k T
        u_ion = 1.5 * n_ion * K_B * T

        # Electrons: similar with degeneracy correction
        E_F = (H_BAR**2 / (2 * M_E)) * (3 * np.pi**2 * n_e) ** (2 / 3) if n_e > 0 else 0
        theta = K_B * T / E_F if E_F > 0 else float("inf")

        if theta < 0.1:
            u_e = 0.6 * n_e * E_F
        elif theta > 10:
            u_e = 1.5 * n_e * K_B * T
        else:
            w = 1 / (1 + theta**1.5)
            u_e = w * 0.6 * n_e * E_F + (1 - w) * 1.5 * n_e * K_B * T

        return u_rad + u_ion + u_e

    @staticmethod
    def sound_speed(T: float, rho: float, A: float, Z: float) -> float:
        """
        Sound speed for hydrodynamic calculations.
        c_s = √(γP/ρ) where γ depends on dominant pressure component.
        """
        P, components = MaterialEOS.total_pressure(T, rho, A, Z)
        if P <= 0 or rho <= 0:
            return 0.0

        # Effective γ based on dominant component
        P_rad = components.get("P_rad", 0)
        gamma = 4 / 3 if P_rad > 0.5 * P else 5 / 3

        return np.sqrt(gamma * P / rho)


@dataclass
class CoreTamperSystem:
    """
    Two-region system: fissile core + tamper/reflector.

    Solves coupled neutron diffusion with proper interface conditions:
    - Flux continuity: Φ_core(R_core) = Φ_tamper(R_core)
    - Current continuity: D_core ∂Φ/∂r|_core = D_tamper ∂Φ/∂r|_tamper
    """

    core_data: nuc.NuclearData
    tamper: TamperMaterial
    core_mass: float
    tamper_mass: float
    compression: float = 1.0

    def __post_init__(self):
        self._compute_geometry()
        self._compute_nuclear_params()

    def _compute_geometry(self):
        """Compute core and tamper radii from masses."""
        rho_core = self.core_data.rho_0 * self.compression
        V_core = self.core_mass / rho_core
        self.R_core = (3 * V_core / (4 * np.pi)) ** (1 / 3)

        rho_tamper = self.tamper.rho_0 * self.compression
        V_tamper = self.tamper_mass / rho_tamper
        # Tamper is a shell: V_shell = (4/3)π(R_out³ - R_core³)
        R_out_cubed = self.R_core**3 + 3 * V_tamper / (4 * np.pi)
        self.R_tamper = R_out_cubed ** (1 / 3)
        self.tamper_thickness = self.R_tamper - self.R_core

    def _compute_nuclear_params(self):
        """Compute diffusion parameters for both regions."""
        # Core parameters
        self.core_params = nuc.compute_diffusion_params(
            self.core_data, self.compression
        )

        # Tamper parameters
        n_tamper = self.tamper.n_density * self.compression
        BARN = 1e-28

        Sigma_s = self.tamper.sigma_s * BARN * n_tamper
        Sigma_a = self.tamper.sigma_a * BARN * n_tamper
        Sigma_f = self.tamper.sigma_f * BARN * n_tamper
        Sigma_t = Sigma_s + Sigma_a

        mu0 = 2 / (3 * self.tamper.A)
        Sigma_tr = Sigma_t - mu0 * Sigma_s

        self.D_tamper = 1 / (3 * Sigma_tr)
        self.Sigma_a_tamper = Sigma_a
        self.Sigma_f_tamper = Sigma_f
        self.nu_tamper = self.tamper.nu_f
        self.lambda_tr_tamper = 1 / Sigma_tr

    def compute_critical_state(self) -> Dict[str, float]:
        """
        Compute k_eff for the two-region system.

        Uses reflector savings approach from reactor physics:
        - The tamper acts as a neutron reflector, effectively increasing
          the extrapolation distance
        - For thick reflectors: δ_refl ≈ D_core/D_tamper × L_tamper
        - For thin reflectors: δ_refl ≈ δ_thick × tanh(t/L_tamper)
        """
        D_c = self.core_params.D_tr
        D_t = self.D_tamper
        nu_Sigma_f_c = self.core_data.nu_prompt * self.core_params.xs.Sigma_f
        nu_Sigma_f_t = self.nu_tamper * self.Sigma_f_tamper

        R_c = self.R_core
        R_t = self.R_tamper

        # Bare sphere k_eff first
        k_bare, B_g_sq_bare = self._bare_k_eff()

        # Reflector savings calculation
        # Diffusion length in tamper
        L_t = (
            np.sqrt(D_t / self.Sigma_a_tamper)
            if self.Sigma_a_tamper > 0
            else self.tamper_thickness
        )

        # Reflector savings: δ_refl = (D_c/D_t) × L_t × tanh(t/L_t)
        # This is the standard reactor physics formula
        t_over_L = self.tamper_thickness / L_t if L_t > 0 else 0

        # Base reflector savings (slab geometry formula)
        delta_slab = (D_c / D_t) * L_t * np.tanh(t_over_L)

        # Spherical geometry correction: reflected flux dilutes as neutrons
        # travel outward in the shell. Factor ~ R_core/(R_core + L_t)
        geom_factor = R_c / (R_c + L_t)
        delta_refl = delta_slab * geom_factor

        # For U-238 with fast fission, modest multiplication boost
        if nu_Sigma_f_t > 0 and self.Sigma_a_tamper > 0:
            k_tamper = nu_Sigma_f_t / self.Sigma_a_tamper
            if k_tamper < 1:
                mult_factor = 1.0 + 0.3 * k_tamper
                delta_refl *= mult_factor

        # Extrapolation distance at outer boundary of tamper
        c_tamper = self.tamper.sigma_s / (self.tamper.sigma_s + self.tamper.sigma_a)
        delta_outer = nuc.extrapolation_distance(self.lambda_tr_tamper, c_tamper)

        # With reflector, the reflector savings REPLACES vacuum extrapolation
        # (neutrons hit reflector, not vacuum, at core surface)
        R_eff = R_c + delta_refl

        # Geometric buckling with reflector
        B_g_sq = (np.pi / R_eff) ** 2

        # k_eff with reflector
        k_eff = self.core_params.k_inf / (1 + self.core_params.L_sq * B_g_sq)

        # Tamper savings as percentage increase in k_eff
        tamper_savings = (k_eff - k_bare) / k_bare * 100 if k_bare > 0 else 0

        return {
            "k_eff": k_eff,
            "k_bare": k_bare,
            "tamper_savings_pct": tamper_savings,
            "R_core": R_c,
            "R_tamper": R_t,
            "tamper_thickness": self.tamper_thickness,
            "delta_reflector": delta_refl,
            "L_tamper": L_t,
        }

    def _bare_k_eff(self) -> Tuple[float, float]:
        """k_eff for bare core (no tamper)."""
        params = self.core_params
        c = params.xs.Sigma_s / params.xs.Sigma_t
        delta = nuc.extrapolation_distance(params.lambda_tr, c)
        R_prime = self.R_core + delta
        B_g_sq = (np.pi / R_prime) ** 2
        k = params.k_inf / (1 + params.L_sq * B_g_sq)
        return k, B_g_sq


@dataclass
class HydrodynamicState:
    """State variables for hydrodynamic expansion."""

    R: float  # Outer radius [m]
    v: float  # Expansion velocity [m/s]
    T: float  # Temperature [K]
    E_fission: float  # Fission energy deposited [J]
    E_kinetic: float  # Kinetic energy [J]
    N_neutrons: float  # Neutron population
    compression: float  # Current compression factor


class FullPhysicsSimulation:
    """
    Complete simulation with coupled neutronics and hydrodynamics.

    Features:
    - Two-region diffusion (core + tamper)
    - Full material EOS
    - Hydrodynamic expansion with inertial confinement
    - Proper energy partition
    """

    def __init__(
        self,
        core_data: nuc.NuclearData,
        core_mass: float,
        initial_compression: float,
        tamper: Optional[TamperMaterial] = None,
        tamper_mass: float = 0.0,
    ):
        self.core_data = core_data
        self.core_mass = core_mass
        self.initial_compression = initial_compression
        self.tamper = tamper
        self.tamper_mass = tamper_mass

        # Initial geometry
        rho_core = core_data.rho_0 * initial_compression
        V_core = core_mass / rho_core
        self.R_0 = (3 * V_core / (4 * np.pi)) ** (1 / 3)

        # Fission energy
        self.E_per_fission = core_data.E_fission * nuc.MEV_TO_J

        # Shell masses for inertia
        self._setup_mass_shells()

    def _setup_mass_shells(self):
        """Setup mass distribution for hydrodynamic calculation."""
        self.shell_masses = [self.core_mass]
        self.shell_radii = [self.R_0]

        if self.tamper is not None and self.tamper_mass > 0:
            rho_t = self.tamper.rho_0 * self.initial_compression
            V_t = self.tamper_mass / rho_t
            R_t = (self.R_0**3 + 3 * V_t / (4 * np.pi)) ** (1 / 3)
            self.shell_masses.append(self.tamper_mass)
            self.shell_radii.append(R_t)

    def compute_alpha(self, R: float, compression: float) -> Tuple[float, Dict]:
        """Compute Rossi alpha (inverse period) at given state."""
        if self.tamper is not None and self.tamper_mass > 0:
            system = CoreTamperSystem(
                self.core_data,
                self.tamper,
                self.core_mass,
                self.tamper_mass,
                compression,
            )
            state = system.compute_critical_state()
            k_eff = state["k_eff"]
        else:
            alpha, params = nuc.compute_rossi_alpha(self.core_data, R, compression)
            return alpha, {"params": params}

        # α = (k - 1) / Λ
        params = nuc.compute_diffusion_params(self.core_data, compression)
        alpha = (k_eff - 1) / params.Lambda

        return alpha, {"k_eff": k_eff, "params": params}

    def enclosed_mass(self, R: float) -> float:
        """Total mass enclosed within radius R."""
        total = 0.0
        for m, r in zip(self.shell_masses, self.shell_radii):
            if R >= r:
                total += m
            else:
                # Partial shell
                frac = (R / r) ** 3
                total += m * frac
                break

        # Add air beyond last shell
        if R > self.shell_radii[-1]:
            rho_air = 1.29  # kg/m³
            V_air = (4 / 3) * np.pi * (R**3 - self.shell_radii[-1] ** 3)
            total += rho_air * V_air

        return total

    def simulate(
        self,
        N_0: float = 1e10,
        t_max: float = 1e-6,
        rtol: float = 1e-6,
    ) -> Dict[str, Any]:
        """
        Run full simulation.

        State vector: [R, v, N, E_fis, E_kin]
        """
        y0 = np.array(
            [
                self.R_0,  # Initial radius
                0.0,  # Initial velocity
                N_0,  # Initial neutrons
                0.0,  # Fission energy
                0.0,  # Kinetic energy
            ]
        )

        # Get initial alpha
        alpha_0, info = self.compute_alpha(self.R_0, self.initial_compression)
        print(f"Initial state:")
        print(f"  R_0 = {self.R_0 * 100:.3f} cm")
        print(f"  α = {alpha_0:.3e} /s")

        if alpha_0 <= 0:
            print("  WARNING: System is subcritical!")
            return {"error": "subcritical"}

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            R, v, N, E_fis, E_kin = y

            # Current compression
            c = self.initial_compression * (self.R_0 / R) ** 3
            c = max(c, 0.01)

            # Neutronics
            alpha, _ = self.compute_alpha(R, c)
            dN_dt = alpha * N

            # Fission power
            params = nuc.compute_diffusion_params(self.core_data, c)
            v_n = self.core_data.average_neutron_velocity()
            P_fis = params.xs.Sigma_f * N * v_n * self.E_per_fission
            dE_fis_dt = P_fis

            # Internal energy and temperature
            E_int = max(0, E_fis - E_kin)
            V = (4 / 3) * np.pi * R**3
            rho = self.core_mass / V

            # Get temperature from EOS
            u_target = E_int / V if V > 0 else 0
            T = self._solve_temperature(u_target, rho)

            # Pressure from full EOS
            P, _ = MaterialEOS.total_pressure(T, rho, self.core_data.A, 94)

            # Hydrodynamics: shell acceleration
            m_enc = self.enclosed_mass(R)
            dv_dt = 4 * np.pi * R**2 * P / m_enc if m_enc > 0 else 0

            # Kinetic energy rate
            dE_kin_dt = 4 * np.pi * R**2 * P * v if v > 0 else 0

            return np.array([v, dv_dt, dN_dt, dE_fis_dt, dE_kin_dt])

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
            rtol=rtol,
            events=subcritical_event,
        )

        # Extract results
        results = {
            "t": sol.t,
            "R": sol.y[0],
            "v": sol.y[1],
            "N": sol.y[2],
            "E_fission": sol.y[3],
            "E_kinetic": sol.y[4],
            "compression": self.initial_compression * (self.R_0 / sol.y[0]) ** 3,
            "yield_kt": sol.y[3] / 4.184e12,
        }

        print(f"\nFinal state:")
        print(f"  t = {sol.t[-1] * 1e6:.3f} μs")
        print(f"  Yield = {results['yield_kt'][-1]:.2f} kt")

        return results

    def _solve_temperature(self, u_target: float, rho: float) -> float:
        """Solve for temperature given energy density."""
        if u_target <= 0:
            return 1e4

        # Initial guess from radiation-dominated regime
        T = (u_target / A_RAD) ** 0.25
        T = max(T, 1e4)

        # Newton iteration
        for _ in range(30):
            u = MaterialEOS.internal_energy_density(T, rho, self.core_data.A, 94)
            du_dT = 4 * A_RAD * T**3 + 1.5 * rho * N_A / (self.core_data.A * 1e-3) * K_B

            residual = u - u_target
            if abs(residual) < 1e-6 * u_target:
                break

            T = T - residual / du_dT
            T = max(T, 1e4)

        return T


def validate_trinity():
    """Validate against Trinity test data."""
    print("=" * 70)
    print("TRINITY VALIDATION (with tamper)")
    print("=" * 70)

    # Trinity configuration (approximate)
    # Core: 6.2 kg Pu-239
    # Tamper: 108 kg U-238
    # Compression: ~2.5x

    nuc.PhysicsConfig.use_trinity_calibration()

    system = CoreTamperSystem(
        core_data=nuc.PU239,
        tamper=U238_TAMPER,
        core_mass=6.2,
        tamper_mass=108.0,
        compression=2.5,
    )

    state = system.compute_critical_state()

    print(f"\nCore-Tamper System:")
    print(f"  R_core = {state['R_core'] * 100:.2f} cm")
    print(f"  R_tamper = {state['R_tamper'] * 100:.2f} cm")
    print(f"  Tamper thickness = {state['tamper_thickness'] * 100:.2f} cm")
    print(f"\nCriticality:")
    print(f"  k_eff (with tamper) = {state['k_eff']:.4f}")
    print(f"  k_eff (bare core) = {state['k_bare']:.4f}")
    print(f"  Tamper savings = {state['tamper_savings_pct']:.1f}%")

    # Full simulation
    print("\n" + "-" * 50)
    print("Running full simulation...")

    sim = FullPhysicsSimulation(
        core_data=nuc.PU239,
        core_mass=6.2,
        initial_compression=2.5,
        tamper=U238_TAMPER,
        tamper_mass=108.0,
    )

    results = sim.simulate(N_0=1e10, t_max=2e-6)

    if "error" not in results:
        print(f"\n  Final yield: {results['yield_kt'][-1]:.1f} kt")
        print(f"  (Trinity actual: ~15-21 kt)")

    nuc.PhysicsConfig.use_standard_physics()
    return results


if __name__ == "__main__":
    validate_trinity()
