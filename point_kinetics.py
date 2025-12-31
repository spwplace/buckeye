"""
Point Kinetics Model - arXiv:1606.01670v1 (Aste 2016)
Neutron diffusion and multiplicative process simulation.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import N_A, c as speed_of_light, Stefan_Boltzmann
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple, Dict, Any, Optional

BARN_TO_M2 = 1e-28
NEUTRON_VELOCITY = 18e6
ENERGY_PER_FISSION_J = 200e6 * 1.602e-19


@dataclass
class NuclearData:
    """Fission spectrum-averaged cross-sections (JEFF database)."""

    name: str
    atomic_mass: float
    isotope_mass_u: float
    density_kg_m3: float
    sigma_f: float
    sigma_s: float
    sigma_c: float
    nu_prompt: float

    @property
    def sigma_t(self) -> float:
        return self.sigma_f + self.sigma_s + self.sigma_c

    @property
    def sigma_a(self) -> float:
        return self.sigma_f + self.sigma_c

    @property
    def atom_density(self) -> float:
        return self.density_kg_m3 / (self.isotope_mass_u * 1.66054e-27)

    @property
    def mu0(self) -> float:
        """Average cosine of scattering angle: mu0 = 2/(3A)"""
        return 2.0 / (3.0 * self.atomic_mass)


PU239 = NuclearData(
    name="Pu-239",
    atomic_mass=239.0,
    isotope_mass_u=239.0521634,
    density_kg_m3=19.86e3,
    sigma_f=1.800,
    sigma_s=4.566 + 1.369,
    sigma_c=0.065,
    nu_prompt=3.165,
)

U235 = NuclearData(
    name="U-235",
    atomic_mass=235.0,
    isotope_mass_u=235.0439299,
    density_kg_m3=18.9e3,
    sigma_f=1.235,
    sigma_s=4.409 + 1.917,
    sigma_c=0.089,
    nu_prompt=2.620,
)


def compute_macroscopic_cross_sections(
    data: NuclearData, compression: float = 1.0
) -> Dict[str, float]:
    """Sigma = sigma * n_atoms, where n_atoms scales with compression."""
    n_atoms = data.atom_density * compression

    return {
        "Sigma_f": data.sigma_f * BARN_TO_M2 * n_atoms,
        "Sigma_s": data.sigma_s * BARN_TO_M2 * n_atoms,
        "Sigma_a": data.sigma_a * BARN_TO_M2 * n_atoms,
        "Sigma_t": data.sigma_t * BARN_TO_M2 * n_atoms,
        "n_atoms": n_atoms,
    }


def compute_diffusion_params(
    data: NuclearData, compression: float = 1.0
) -> Dict[str, float]:
    """Transport MFP, diffusion coefficient, material buckling."""
    xs = compute_macroscopic_cross_sections(data, compression)

    # Sigma_tr = Sigma_t - mu0 * Sigma_s (Eq. 22)
    Sigma_tr = xs["Sigma_t"] - data.mu0 * xs["Sigma_s"]
    lambda_tr = 1.0 / Sigma_tr

    # D = lambda_tr / 3
    D = lambda_tr / 3.0

    # B_m^2 = (nu*Sigma_f - Sigma_a) / D (Eq. 14)
    B_m_squared = (data.nu_prompt * xs["Sigma_f"] - xs["Sigma_a"]) / D

    L_squared = D / xs["Sigma_a"]
    k_inf = data.nu_prompt * xs["Sigma_f"] / xs["Sigma_a"]

    return {
        "lambda_tr": lambda_tr,
        "D": D,
        "B_m_squared": B_m_squared,
        "L_squared": L_squared,
        "k_inf": k_inf,
        **xs,
    }


def compute_critical_radius(
    data: NuclearData, compression: float = 1.0
) -> Tuple[Optional[float], Optional[float], Dict[str, float]]:
    """
    Critical radius via diffusion theory.
    Flux vanishes at R' = R + 0.71045*lambda_tr (Eq. 20)
    At criticality: B_g^2 = B_m^2 where B_g^2 = pi^2/R'^2 (Eq. 21)
    """
    params = compute_diffusion_params(data, compression)
    delta = 0.71045 * params["lambda_tr"]

    if params["B_m_squared"] <= 0:
        return None, None, params

    R_prime_c = np.pi / np.sqrt(params["B_m_squared"])
    R_c = R_prime_c - delta

    return R_c, delta, params


def compute_critical_mass(
    data: NuclearData, compression: float = 1.0
) -> Tuple[Optional[float], Optional[float], Dict[str, float]]:
    """Critical mass for bare sphere: m = (4/3)*pi*R^3*rho*c"""
    R_c, delta, params = compute_critical_radius(data, compression)

    if R_c is None or R_c <= 0:
        return None, None, params

    volume = (4.0 / 3.0) * np.pi * R_c**3
    mass = volume * data.density_kg_m3 * compression

    return mass, R_c, params


@dataclass
class CoreState:
    radius: float
    velocity: float
    neutron_count: float
    fission_energy: float
    kinetic_energy: float


def inverse_bomb_period(
    data: NuclearData,
    radius: float,
    compression: float,
    extrapolation_factor: float = 0.71045,
) -> Tuple[float, Dict[str, float]]:
    """
    Rossi alpha: alpha = v*D*(B_m^2 - B_g^2) (Eq. 13)
    alpha > 0: supercritical, alpha < 0: subcritical
    """
    params = compute_diffusion_params(data, compression)

    delta = extrapolation_factor * params["lambda_tr"]
    R_prime = radius + delta

    # B_g^2 = pi^2 / R'^2 (Eq. 16)
    B_g_squared = (np.pi / R_prime) ** 2

    alpha = NEUTRON_VELOCITY * params["D"] * (params["B_m_squared"] - B_g_squared)

    return alpha, params


def mass_profile(
    radius: float, R_min: float, compression: float, core_mass: float = 6.2
) -> float:
    """Enclosed mass vs radius (simplified: core + air)."""
    if radius <= R_min:
        return core_mass

    R_0 = R_min * compression ** (1 / 3)
    if radius > R_0:
        air_density = 1.29
        air_volume = (4 / 3) * np.pi * (radius**3 - R_0**3)
        return core_mass + air_density * air_volume

    return core_mass


def point_kinetics_rhs(
    t: float,
    state: np.ndarray,
    data: NuclearData,
    initial_radius: float,
    compression: float,
    core_mass: float,
) -> np.ndarray:
    """ODE RHS: state = [R, dR/dt, N, E_fission, E_kinetic]"""
    R, R_dot, N, E_fission, E_kinetic = state

    current_compression = (initial_radius / R) ** 3 * compression
    if current_compression < 0.01:
        current_compression = 0.01

    alpha, params = inverse_bomb_period(data, R, current_compression)

    # dN/dt = alpha * N (Eq. 17)
    dN_dt = alpha * N

    # P = Sigma_f * N * v * epsilon_f (Eq. 31)
    power = params["Sigma_f"] * N * NEUTRON_VELOCITY * ENERGY_PER_FISSION_J
    dE_fission_dt = power

    E_radiation = max(0.0, E_fission - E_kinetic)

    # p = E_rad / (3V) (Eq. 33)
    volume = (4 / 3) * np.pi * R**3
    pressure = E_radiation / (3 * volume) if volume > 0 else 0

    # d(R_dot)/dt = 4*pi*R^2 * p / m (from Eq. 35)
    m_eff = mass_profile(R, initial_radius, compression, core_mass)
    if R_dot > 0 and m_eff > 0:
        dR_dot_dt = 4 * np.pi * R**2 * pressure / m_eff
    else:
        dR_dot_dt = 4 * np.pi * R**2 * pressure / core_mass if pressure > 0 else 0

    dE_kinetic_dt = 4 * np.pi * R**2 * pressure * R_dot if R_dot > 0 else 0

    return np.array([R_dot, dR_dot_dt, dN_dt, dE_fission_dt, dE_kinetic_dt])


def simulate_explosion(
    data: NuclearData,
    core_mass_kg: float,
    compression: float,
    initial_neutrons: float = 1e8,
    t_max: float = 1e-6,
    max_step: float = 1e-12,
) -> Dict[str, Any]:
    """Simulate prompt supercritical excursion via point kinetics."""

    compressed_density = data.density_kg_m3 * compression
    initial_volume = core_mass_kg / compressed_density
    initial_radius = (3 * initial_volume / (4 * np.pi)) ** (1 / 3)

    alpha_0, params_0 = inverse_bomb_period(data, initial_radius, compression)
    k_eff_0 = 1 + alpha_0 * (
        1 / (data.nu_prompt * params_0["Sigma_f"] * NEUTRON_VELOCITY)
    )

    print(f"Initial conditions:")
    print(f"  Radius: {initial_radius * 100:.3f} cm")
    print(f"  Compression: {compression:.2f}x")
    print(f"  Alpha: {alpha_0:.3e} /s")
    print(f"  k_eff: {k_eff_0:.4f}")
    print(f"  Supercritical: {alpha_0 > 0}")

    y0 = np.array([initial_radius, 0.0, initial_neutrons, 0.0, 0.0])

    def subcritical_event(t: float, y: np.ndarray) -> float:
        R = y[0]
        current_compression = (initial_radius / R) ** 3 * compression
        if current_compression < 0.01:
            return -1.0
        alpha, _ = inverse_bomb_period(data, R, current_compression)
        return alpha

    subcritical_event.terminal = True  # type: ignore[attr-defined]
    subcritical_event.direction = -1  # type: ignore[attr-defined]

    solution = solve_ivp(
        lambda t, y: point_kinetics_rhs(
            t, y, data, initial_radius, compression, core_mass_kg
        ),
        (0, t_max),
        y0,
        method="RK45",
        max_step=max_step,
        events=subcritical_event,
        dense_output=True,
    )

    t = solution.t
    R = solution.y[0]
    R_dot = solution.y[1]
    N = solution.y[2]
    E_fission = solution.y[3]
    E_kinetic = solution.y[4]

    compression_history = (initial_radius / R) ** 3 * compression
    E_radiation = np.maximum(0, E_fission - E_kinetic)

    volume = (4 / 3) * np.pi * R**3
    pressure = E_radiation / (3 * volume)

    # T = (3pc / 4*sigma)^(1/4) from radiation pressure relation
    T_radiation = (3 * pressure * speed_of_light / (4 * Stefan_Boltzmann)) ** 0.25

    alpha_history = np.zeros_like(t)
    for i in range(len(t)):
        alpha_history[i], _ = inverse_bomb_period(data, R[i], compression_history[i])

    yield_kt = E_fission / 4.184e12

    theoretical_max = (
        core_mass_kg / (data.isotope_mass_u * 1.66054e-27) * ENERGY_PER_FISSION_J
    )
    efficiency = E_fission / theoretical_max

    results = {
        "t": t,
        "t_shakes": t * 1e8,
        "R": R,
        "R_cm": R * 100,
        "R_dot": R_dot,
        "R_dot_km_s": R_dot / 1000,
        "N": N,
        "E_fission_J": E_fission,
        "E_kinetic_J": E_kinetic,
        "E_radiation_J": E_radiation,
        "compression": compression_history,
        "pressure_Pa": pressure,
        "pressure_Gbar": pressure / 1e14,
        "T_radiation_K": T_radiation,
        "alpha": alpha_history,
        "yield_kt": yield_kt,
        "efficiency": efficiency,
        "final_yield_kt": yield_kt[-1] if len(yield_kt) > 0 else 0,
        "final_efficiency": efficiency[-1] if len(efficiency) > 0 else 0,
        "solution": solution,
    }

    print(f"\nResults:")
    print(f"  Final yield: {results['final_yield_kt']:.2f} kt")
    print(f"  Efficiency: {results['final_efficiency'] * 100:.2f}%")
    print(f"  Peak temperature: {np.max(T_radiation):.2e} K")
    print(f"  Peak pressure: {np.max(pressure) / 1e14:.2f} Gbar")
    print(f"  Time: {t[-1] * 1e6:.3f} us")

    return results


def plot_results(results: Dict[str, Any], title: str = "Point Kinetics"):
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle(title, fontsize=14)

    t_us = results["t"] * 1e6

    axes[0, 0].semilogy(t_us, results["N"])
    axes[0, 0].set_xlabel("Time (us)")
    axes[0, 0].set_ylabel("Neutron count")
    axes[0, 0].set_title("Neutron Population")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(t_us, results["yield_kt"])
    axes[0, 1].set_xlabel("Time (us)")
    axes[0, 1].set_ylabel("Energy (kt)")
    axes[0, 1].set_title("Energy Release")
    axes[0, 1].grid(True, alpha=0.3)

    axes[0, 2].plot(t_us, results["alpha"] / 1e8, "b-")
    axes[0, 2].axhline(y=0, color="r", linestyle="--", alpha=0.5)
    axes[0, 2].set_xlabel("Time (us)")
    axes[0, 2].set_ylabel("Alpha (10^8 /s)")
    axes[0, 2].set_title("Rossi Alpha")
    axes[0, 2].grid(True, alpha=0.3)

    axes[1, 0].plot(t_us, results["R_cm"], "b-")
    axes[1, 0].set_xlabel("Time (us)")
    axes[1, 0].set_ylabel("Radius (cm)", color="b")
    ax2 = axes[1, 0].twinx()
    ax2.plot(t_us, results["R_dot_km_s"], "r-")
    ax2.set_ylabel("Velocity (km/s)", color="r")
    axes[1, 0].set_title("Core Expansion")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(t_us, results["T_radiation_K"])
    axes[1, 1].set_xlabel("Time (us)")
    axes[1, 1].set_ylabel("Temperature (K)")
    axes[1, 1].set_title("Radiation Temperature")
    axes[1, 1].grid(True, alpha=0.3)

    axes[1, 2].plot(t_us, results["pressure_Gbar"])
    axes[1, 2].set_xlabel("Time (us)")
    axes[1, 2].set_ylabel("Pressure (Gbar)")
    axes[1, 2].set_title("Radiation Pressure")
    axes[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    print("=" * 60)
    print("Critical Mass Calculations")
    print("=" * 60)

    for material in [PU239, U235]:
        print(f"\n{material.name}:")
        mass, R_c, params = compute_critical_mass(material)
        if mass is not None and R_c is not None:
            print(f"  R_c: {R_c * 100:.3f} cm")
            print(f"  m_c: {mass:.2f} kg")
            print(f"  k_inf: {params['k_inf']:.4f}")

        print(f"\n  Compression effects:")
        for c in [1.0, 1.5, 2.0, 2.5]:
            mass_c, R_c_c, _ = compute_critical_mass(material, compression=c)
            if mass_c is not None and R_c_c is not None:
                print(f"    c={c:.1f}: R_c={R_c_c * 100:.2f} cm, m_c={mass_c:.2f} kg")

    print("\n" + "=" * 60)
    print("Simulation")
    print("=" * 60)

    results = simulate_explosion(
        data=PU239,
        core_mass_kg=6.2,
        compression=2.5,
        initial_neutrons=1e8,
        t_max=2e-6,
        max_step=1e-12,
    )

    fig = plot_results(results, "Point Kinetics - Pu-239 (c=2.5)")
    plt.savefig("point_kinetics_results.png", dpi=150, bbox_inches="tight")
    plt.close()

    print("\n" + "=" * 60)
    print("Compression Study")
    print("=" * 60)

    compressions = [2.2, 2.4, 2.6, 2.8]
    yields = []

    for c in compressions:
        print(f"\nc={c}:")
        res = simulate_explosion(
            data=PU239,
            core_mass_kg=6.2,
            compression=c,
            initial_neutrons=1e8,
            t_max=2e-6,
            max_step=1e-12,
        )
        yields.append(res["final_yield_kt"])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(compressions, yields, "bo-", markersize=8)
    ax.set_xlabel("Compression Factor")
    ax.set_ylabel("Yield (kt)")
    ax.set_title("Yield vs Compression")
    ax.grid(True, alpha=0.3)
    plt.savefig("yield_vs_compression.png", dpi=150, bbox_inches="tight")
    plt.close()

    print("\nPlots saved: point_kinetics_results.png, yield_vs_compression.png")
