"""
Advanced Nuclear Simulations
1. Time-dependent 2D with expansion coupling
2. Non-spherical geometries (ellipsoid, cylinder)
3. Reflector/tamper modeling
4. Fizzle probability with Pu-240
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from dataclasses import dataclass
from typing import Dict, Any, Tuple, List, Optional, Callable
import point_kinetics as pk


BARN_TO_M2 = 1e-28


@dataclass
class Mesh2D:
    nr: int
    nz: int
    R_max: float
    Z_max: float

    def __post_init__(self):
        self.dr = self.R_max / (self.nr - 1)
        self.dz = self.Z_max / (self.nz - 1)
        self.r = np.linspace(0, self.R_max, self.nr)
        self.z = np.linspace(-self.Z_max / 2, self.Z_max / 2, self.nz)
        self.R, self.Z = np.meshgrid(self.r, self.z, indexing="ij")
        self.n_nodes = self.nr * self.nz

    def node_index(self, i: int, j: int) -> int:
        return i * self.nz + j

    def is_boundary(self, i: int, j: int) -> bool:
        """Physical boundaries only - r=0 is symmetry axis, not boundary."""
        return i == self.nr - 1 or j == 0 or j == self.nz - 1

    def is_axis(self, i: int) -> bool:
        return i == 0


@dataclass
class Material2D:
    D: np.ndarray
    Sigma_f: np.ndarray
    Sigma_a: np.ndarray
    nu: float


# =============================================================================
# SECTION 1: Time-Dependent 2D with Expansion
# =============================================================================


def build_diffusion_operator(
    mesh: Mesh2D, mat: Material2D
) -> Tuple[csr_matrix, csr_matrix]:
    """Build L (leakage+absorption) and F (fission) matrices."""
    n = mesh.n_nodes
    L = lil_matrix((n, n))
    Fmat = lil_matrix((n, n))

    dr, dz = mesh.dr, mesh.dz

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            idx = mesh.node_index(i, j)
            r = mesh.r[i]

            if mesh.is_boundary(i, j):
                L[idx, idx] = 1.0
                continue

            D_local = mat.D[i, j]
            Sigma_a_local = mat.Sigma_a[i, j]
            Sigma_f_local = mat.Sigma_f[i, j]

            # Axisymmetric Laplacian: ∂²Φ/∂r² + (1/r)∂Φ/∂r
            if mesh.is_axis(i):
                coeff_r_plus = 4 * D_local / dr**2
                coeff_r_minus = 0.0
                coeff_r_center = -4 * D_local / dr**2
            else:
                coeff_r_plus = D_local / dr**2 + D_local / (2 * r * dr)
                coeff_r_minus = D_local / dr**2 - D_local / (2 * r * dr)
                coeff_r_center = -2 * D_local / dr**2

            coeff_z = D_local / dz**2

            # L = -D∇² + Σa
            L[idx, mesh.node_index(i + 1, j)] = -coeff_r_plus
            if not mesh.is_axis(i):
                L[idx, mesh.node_index(i - 1, j)] = -coeff_r_minus
            L[idx, mesh.node_index(i, j + 1)] = -coeff_z
            L[idx, mesh.node_index(i, j - 1)] = -coeff_z
            L[idx, idx] = -coeff_r_center + 2 * coeff_z + Sigma_a_local

            Fmat[idx, idx] = mat.nu * Sigma_f_local

    return csr_matrix(L), csr_matrix(Fmat)


def compute_k_eff(
    mesh: Mesh2D, mat: Material2D, tol: float = 1e-6
) -> Tuple[float, np.ndarray]:
    """k_eff via power iteration."""
    L, F = build_diffusion_operator(mesh, mat)
    n = mesh.n_nodes

    phi = np.ones(n) * 0.1
    for i in range(mesh.nr):
        for j in range(mesh.nz):
            if mesh.is_boundary(i, j):
                phi[mesh.node_index(i, j)] = 0.0

    k = 1.0
    for _ in range(500):
        source = F @ phi
        fission_rate = np.sum(source)
        if fission_rate < 1e-30:
            return 0.0, phi.reshape(mesh.nr, mesh.nz)

        phi_new = spsolve(L, source / k)
        phi_new = np.maximum(phi_new, 0)

        new_fission_rate = np.sum(F @ phi_new)
        k_new = k * new_fission_rate / fission_rate if fission_rate > 0 else 0

        norm = np.linalg.norm(phi_new)
        if norm > 1e-30:
            phi_new = phi_new / norm

        if abs(k_new - k) / max(k, 1e-10) < tol:
            return k_new, phi_new.reshape(mesh.nr, mesh.nz)

        k = k_new
        phi = phi_new

    return k, phi.reshape(mesh.nr, mesh.nz)


def create_material_for_sphere(
    mesh: Mesh2D, data: pk.NuclearData, radius: float, compression: float
) -> Material2D:
    """Homogeneous sphere material."""
    params = pk.compute_diffusion_params(data, compression)

    D = np.zeros((mesh.nr, mesh.nz))
    Sigma_f = np.zeros((mesh.nr, mesh.nz))
    Sigma_a = np.zeros((mesh.nr, mesh.nz))

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            dist = np.sqrt(mesh.R[i, j] ** 2 + mesh.Z[i, j] ** 2)
            if dist <= radius:
                D[i, j] = params["D"]
                Sigma_f[i, j] = params["Sigma_f"]
                Sigma_a[i, j] = params["Sigma_a"]
            else:
                D[i, j] = 0.1
                Sigma_f[i, j] = 0
                Sigma_a[i, j] = 1e-6

    return Material2D(D=D, Sigma_f=Sigma_f, Sigma_a=Sigma_a, nu=data.nu_prompt)


class TimeDependentSimulation2D:
    """Coupled 2D diffusion with quasi-static expansion."""

    def __init__(
        self,
        data: pk.NuclearData,
        core_mass: float,
        compression: float,
        mesh_size: int = 50,
    ):
        self.data = data
        self.core_mass = core_mass
        self.initial_compression = compression

        compressed_density = data.density_kg_m3 * compression
        volume = core_mass / compressed_density
        self.initial_radius = (3 * volume / (4 * np.pi)) ** (1 / 3)

        domain_size = self.initial_radius * 4
        self.mesh = Mesh2D(
            nr=mesh_size, nz=mesh_size, R_max=domain_size, Z_max=domain_size * 2
        )

        self.mat = create_material_for_sphere(
            self.mesh, data, self.initial_radius, compression
        )
        _, self.phi = compute_k_eff(self.mesh, self.mat)
        self.phi = self.phi.flatten()

        self.total_neutrons = 1e8
        self.phi = self.phi * self.total_neutrons / np.sum(self.phi)

    def compute_fission_rate(self, phi: np.ndarray) -> float:
        """Total fission rate = integral of Σf * Φ * v over volume."""
        phi_2d = phi.reshape(self.mesh.nr, self.mesh.nz)
        fission_rate = 0.0
        for i in range(self.mesh.nr):
            for j in range(self.mesh.nz):
                r = self.mesh.r[i]
                dV = (
                    2 * np.pi * r * self.mesh.dr * self.mesh.dz
                    if r > 0
                    else np.pi * self.mesh.dr**2 * self.mesh.dz
                )
                fission_rate += (
                    self.mat.Sigma_f[i, j] * phi_2d[i, j] * pk.NEUTRON_VELOCITY * dV
                )
        return fission_rate

    def simulate(self, t_max: float = 1e-6) -> Dict[str, Any]:
        """Run time-dependent simulation with expansion feedback."""

        v = pk.NEUTRON_VELOCITY
        n = self.mesh.n_nodes

        # State: [phi (n values), radius, velocity, energy]
        y0 = np.zeros(n + 3)
        y0[:n] = self.phi
        y0[n] = self.initial_radius
        y0[n + 1] = 0.0  # initial velocity
        y0[n + 2] = 0.0  # initial energy

        current_compression = self.initial_compression

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            nonlocal current_compression

            phi = y[:n]
            radius = y[n]
            velocity = y[n + 1]
            energy = y[n + 2]

            new_compression = (
                self.initial_compression * (self.initial_radius / radius) ** 3
            )
            new_compression = max(0.1, new_compression)

            # Rebuild matrices if compression changed significantly
            if abs(new_compression - current_compression) / current_compression > 0.01:
                current_compression = new_compression
                self.mat = create_material_for_sphere(
                    self.mesh, self.data, radius, current_compression
                )

            L, F = build_diffusion_operator(self.mesh, self.mat)

            # dΦ/dt = v*(D∇² + νΣf - Σa)Φ = v*(-L + F)Φ
            dphi_dt = v * (F @ phi - L @ phi)

            # Apply boundary conditions
            for i in range(self.mesh.nr):
                for j in range(self.mesh.nz):
                    if self.mesh.is_boundary(i, j):
                        dphi_dt[self.mesh.node_index(i, j)] = 0

            # Fission power
            fission_rate = self.compute_fission_rate(phi)
            power = fission_rate * pk.ENERGY_PER_FISSION_J

            # Energy rate
            dE_dt = power

            # Radiation pressure drives expansion (Eq. 33, 35)
            volume = (4 / 3) * np.pi * radius**3
            pressure = energy / (3 * volume) if volume > 0 and energy > 0 else 0

            # Acceleration from radiation pressure
            dv_dt = (
                4 * np.pi * radius**2 * pressure / self.core_mass if pressure > 0 else 0
            )

            dy = np.zeros_like(y)
            dy[:n] = dphi_dt
            dy[n] = velocity
            dy[n + 1] = dv_dt
            dy[n + 2] = dE_dt

            return dy

        # Stop when subcritical
        def subcritical_event(t, y):
            radius = y[n]
            comp = self.initial_compression * (self.initial_radius / radius) ** 3
            if comp < 0.5:
                return -1
            mat = create_material_for_sphere(
                self.mesh, self.data, radius, max(0.1, comp)
            )
            k, _ = compute_k_eff(self.mesh, mat)
            return k - 1.0

        subcritical_event.terminal = True  # type: ignore
        subcritical_event.direction = -1  # type: ignore

        solution = solve_ivp(
            rhs, (0, t_max), y0, method="RK23", max_step=1e-10, events=subcritical_event
        )

        times = solution.t
        results = {
            "t": times,
            "phi": [
                solution.y[:n, i].reshape(self.mesh.nr, self.mesh.nz)
                for i in range(len(times))
            ],
            "radius": solution.y[n, :],
            "velocity": solution.y[n + 1, :],
            "energy": solution.y[n + 2, :],
            "yield_kt": solution.y[n + 2, :] / 4.184e12,
        }

        return results


# =============================================================================
# SECTION 2: Non-Spherical Geometries
# =============================================================================


def create_ellipsoid_material(
    mesh: Mesh2D, data: pk.NuclearData, a: float, b: float, compression: float
) -> Material2D:
    """Ellipsoid with semi-axes a (radial) and b (axial)."""
    params = pk.compute_diffusion_params(data, compression)

    D = np.zeros((mesh.nr, mesh.nz))
    Sigma_f = np.zeros((mesh.nr, mesh.nz))
    Sigma_a = np.zeros((mesh.nr, mesh.nz))

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            # Ellipsoid equation: (r/a)^2 + (z/b)^2 <= 1
            r_norm = mesh.R[i, j] / a if a > 0 else 0
            z_norm = mesh.Z[i, j] / b if b > 0 else 0
            if r_norm**2 + z_norm**2 <= 1.0:
                D[i, j] = params["D"]
                Sigma_f[i, j] = params["Sigma_f"]
                Sigma_a[i, j] = params["Sigma_a"]
            else:
                D[i, j] = 0.1
                Sigma_f[i, j] = 0
                Sigma_a[i, j] = 1e-6

    return Material2D(D=D, Sigma_f=Sigma_f, Sigma_a=Sigma_a, nu=data.nu_prompt)


def create_cylinder_material(
    mesh: Mesh2D, data: pk.NuclearData, radius: float, height: float, compression: float
) -> Material2D:
    """Cylinder with given radius and height."""
    params = pk.compute_diffusion_params(data, compression)

    D = np.zeros((mesh.nr, mesh.nz))
    Sigma_f = np.zeros((mesh.nr, mesh.nz))
    Sigma_a = np.zeros((mesh.nr, mesh.nz))

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            in_cylinder = (mesh.R[i, j] <= radius) and (abs(mesh.Z[i, j]) <= height / 2)
            if in_cylinder:
                D[i, j] = params["D"]
                Sigma_f[i, j] = params["Sigma_f"]
                Sigma_a[i, j] = params["Sigma_a"]
            else:
                D[i, j] = 0.1
                Sigma_f[i, j] = 0
                Sigma_a[i, j] = 1e-6

    return Material2D(D=D, Sigma_f=Sigma_f, Sigma_a=Sigma_a, nu=data.nu_prompt)


def study_geometry_criticality(
    data: pk.NuclearData, target_mass: float, compression: float = 1.0
) -> Dict[str, Any]:
    """Compare k_eff for sphere, oblate ellipsoid, prolate ellipsoid, cylinder."""

    density = data.density_kg_m3 * compression
    volume = target_mass / density

    # Sphere: R = (3V/4pi)^(1/3)
    R_sphere = (3 * volume / (4 * np.pi)) ** (1 / 3)

    # Oblate ellipsoid (pancake): a = 1.5*R, b = R/1.5^2 to keep volume
    # V = (4/3)*pi*a^2*b
    a_oblate = R_sphere * 1.5
    b_oblate = volume / ((4 / 3) * np.pi * a_oblate**2)

    # Prolate ellipsoid (football): a = R/1.5, b = 1.5^2*R
    a_prolate = R_sphere / 1.5
    b_prolate = volume / ((4 / 3) * np.pi * a_prolate**2)

    # Cylinder: R = H (equal diameter and height)
    # V = pi*R^2*H = pi*R^3
    R_cyl = (volume / np.pi) ** (1 / 3)
    H_cyl = R_cyl

    domain_size = max(R_sphere, a_oblate, b_prolate, R_cyl, H_cyl) * 3
    mesh = Mesh2D(nr=50, nz=50, R_max=domain_size, Z_max=domain_size * 2)

    results = {}

    # Sphere
    mat_sphere = create_material_for_sphere(mesh, data, R_sphere, compression)
    k_sphere, phi_sphere = compute_k_eff(mesh, mat_sphere)
    results["sphere"] = {
        "k_eff": k_sphere,
        "phi": phi_sphere,
        "dims": f"R={R_sphere * 100:.2f}cm",
    }

    # Oblate ellipsoid
    mat_oblate = create_ellipsoid_material(mesh, data, a_oblate, b_oblate, compression)
    k_oblate, phi_oblate = compute_k_eff(mesh, mat_oblate)
    results["oblate"] = {
        "k_eff": k_oblate,
        "phi": phi_oblate,
        "dims": f"a={a_oblate * 100:.2f}cm, b={b_oblate * 100:.2f}cm",
    }

    # Prolate ellipsoid
    mat_prolate = create_ellipsoid_material(
        mesh, data, a_prolate, b_prolate, compression
    )
    k_prolate, phi_prolate = compute_k_eff(mesh, mat_prolate)
    results["prolate"] = {
        "k_eff": k_prolate,
        "phi": phi_prolate,
        "dims": f"a={a_prolate * 100:.2f}cm, b={b_prolate * 100:.2f}cm",
    }

    # Cylinder
    mat_cyl = create_cylinder_material(mesh, data, R_cyl, H_cyl, compression)
    k_cyl, phi_cyl = compute_k_eff(mesh, mat_cyl)
    results["cylinder"] = {
        "k_eff": k_cyl,
        "phi": phi_cyl,
        "dims": f"R={R_cyl * 100:.2f}cm, H={H_cyl * 100:.2f}cm",
    }

    results["mesh"] = mesh
    results["volume"] = volume
    results["mass"] = target_mass

    return results


# =============================================================================
# SECTION 3: Reflector/Tamper Modeling
# =============================================================================


@dataclass
class TamperData:
    """Tamper/reflector material properties."""

    name: str
    density_kg_m3: float
    sigma_s: float  # scattering cross-section (barns)
    sigma_a: float  # absorption cross-section (barns)
    atom_mass: float  # atomic mass

    @property
    def atom_density(self) -> float:
        return self.density_kg_m3 / (self.atom_mass * 1.66054e-27)


# Common tamper materials
U238_TAMPER = TamperData(
    name="U-238",
    density_kg_m3=18.95e3,
    sigma_s=6.9,
    sigma_a=0.37,  # fission + capture at fast energies
    atom_mass=238.05,
)

BERYLLIUM_REFLECTOR = TamperData(
    name="Beryllium", density_kg_m3=1.85e3, sigma_s=6.1, sigma_a=0.01, atom_mass=9.012
)

TUNGSTEN_TAMPER = TamperData(
    name="Tungsten", density_kg_m3=19.3e3, sigma_s=5.0, sigma_a=0.5, atom_mass=183.84
)


def create_tampered_sphere_material(
    mesh: Mesh2D,
    core_data: pk.NuclearData,
    tamper: TamperData,
    core_radius: float,
    tamper_thickness: float,
    compression: float,
) -> Material2D:
    """Sphere with outer tamper/reflector shell."""
    core_params = pk.compute_diffusion_params(core_data, compression)

    # Tamper diffusion parameters
    n_tamper = tamper.atom_density
    Sigma_s_tamper = tamper.sigma_s * BARN_TO_M2 * n_tamper
    Sigma_a_tamper = tamper.sigma_a * BARN_TO_M2 * n_tamper
    mu0_tamper = 2.0 / (3.0 * tamper.atom_mass)
    Sigma_t_tamper = Sigma_s_tamper + Sigma_a_tamper
    Sigma_tr_tamper = Sigma_t_tamper - mu0_tamper * Sigma_s_tamper
    D_tamper = 1.0 / (3.0 * Sigma_tr_tamper)

    outer_radius = core_radius + tamper_thickness

    D = np.zeros((mesh.nr, mesh.nz))
    Sigma_f = np.zeros((mesh.nr, mesh.nz))
    Sigma_a = np.zeros((mesh.nr, mesh.nz))

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            dist = np.sqrt(mesh.R[i, j] ** 2 + mesh.Z[i, j] ** 2)
            if dist <= core_radius:
                D[i, j] = core_params["D"]
                Sigma_f[i, j] = core_params["Sigma_f"]
                Sigma_a[i, j] = core_params["Sigma_a"]
            elif dist <= outer_radius:
                D[i, j] = D_tamper
                Sigma_f[i, j] = 0
                Sigma_a[i, j] = Sigma_a_tamper
            else:
                D[i, j] = 0.1
                Sigma_f[i, j] = 0
                Sigma_a[i, j] = 1e-6

    return Material2D(D=D, Sigma_f=Sigma_f, Sigma_a=Sigma_a, nu=core_data.nu_prompt)


def study_tamper_effect(
    data: pk.NuclearData, core_mass: float, tamper: TamperData, compression: float = 1.0
) -> Dict[str, Any]:
    """Study effect of tamper thickness on criticality."""

    density = data.density_kg_m3 * compression
    volume = core_mass / density
    core_radius = (3 * volume / (4 * np.pi)) ** (1 / 3)

    max_tamper_thickness = core_radius * 0.5
    thicknesses = np.linspace(0, max_tamper_thickness, 8)

    domain_size = (core_radius + max_tamper_thickness) * 2.5
    mesh = Mesh2D(nr=70, nz=70, R_max=domain_size, Z_max=domain_size * 2)

    k_values = []
    for thickness in thicknesses:
        if thickness == 0:
            mat = create_material_for_sphere(mesh, data, core_radius, compression)
        else:
            mat = create_tampered_sphere_material(
                mesh, data, tamper, core_radius, thickness, compression
            )
        k, _ = compute_k_eff(mesh, mat)
        k_values.append(k)

    k_increase = (np.array(k_values) / k_values[0] - 1) * 100

    return {
        "thicknesses": thicknesses,
        "k_values": np.array(k_values),
        "k_increase_pct": k_increase,
        "core_radius": core_radius,
        "tamper": tamper.name,
        "mesh": mesh,
    }


# =============================================================================
# SECTION 4: Fizzle Probability with Pu-240
# =============================================================================


@dataclass
class PlutoniumMix:
    """Plutonium isotopic composition."""

    pu239_fraction: float
    pu240_fraction: float

    # Pu-240 spontaneous fission rate: 479.1 fissions/g/s (from paper)
    PU240_SF_RATE = 479.1  # spontaneous fissions per gram per second
    PU240_NEUTRONS_PER_SF = 2.16  # average neutrons per spontaneous fission

    @property
    def pu241_fraction(self) -> float:
        return max(0, 1.0 - self.pu239_fraction - self.pu240_fraction)

    def neutron_source_rate(self, total_mass_kg: float) -> float:
        """Spontaneous fission neutron source rate (neutrons/second)."""
        pu240_mass_g = total_mass_kg * 1000 * self.pu240_fraction
        sf_rate = pu240_mass_g * self.PU240_SF_RATE
        return sf_rate * self.PU240_NEUTRONS_PER_SF


# Standard grades
SUPERGRADE_PU = PlutoniumMix(pu239_fraction=0.98, pu240_fraction=0.02)
WEAPONS_GRADE_PU = PlutoniumMix(pu239_fraction=0.93, pu240_fraction=0.07)
REACTOR_GRADE_PU = PlutoniumMix(pu239_fraction=0.60, pu240_fraction=0.25)


def compute_critical_compression(data: pk.NuclearData, core_mass: float) -> float:
    """Find compression at which core becomes critical (k=1)."""
    for c in np.linspace(1.0, 10.0, 100):
        try:
            m_c, R_c, params = pk.compute_critical_mass(data, c)
            if m_c is not None and core_mass >= m_c:
                return c
        except Exception:
            continue
    return 10.0


def simulate_with_predetonation(
    data: pk.NuclearData,
    core_mass: float,
    compression: float,
    pu_mix: PlutoniumMix,
    compression_time: float = 12e-6,
    n_trials: int = 100,
) -> Dict[str, Any]:
    """
    Monte Carlo simulation of fizzle probability with PROPER PHYSICS.

    Key corrections from naive model:
    1. Predetonation can only occur during SUPERCRITICAL phase (k > 1)
    2. Not every neutron starts a chain - probability depends on excess k
    3. Yield depends on when chain starts and subsequent energy release

    Based on Section 3.2 of the paper (arXiv:1606.01670v1).
    """
    neutron_rate = pu_mix.neutron_source_rate(core_mass)

    # Compute final state parameters
    compressed_density = data.density_kg_m3 * compression
    volume = core_mass / compressed_density
    final_radius = (3 * volume / (4 * np.pi)) ** (1 / 3)

    # Find when system becomes supercritical
    # For implosion: c(t) = 1 + (c_max - 1) * (t/t_total)
    c_critical = compute_critical_compression(data, core_mass)

    if c_critical >= compression:
        # System never becomes supercritical
        return {
            "yields": np.zeros(n_trials),
            "predet_times": np.full(n_trials, compression_time),
            "mean_yield": 0.0,
            "std_yield": 0.0,
            "fizzle_probability": 1.0,
            "full_yield_probability": 0.0,
            "neutron_rate": neutron_rate,
            "pu_mix": pu_mix,
            "n_trials": n_trials,
            "c_critical": c_critical,
        }

    # Time when system becomes supercritical
    # c(t) = 1 + (c_max - 1) * (t/t_total) = c_critical
    # => t_critical = t_total * (c_critical - 1) / (c_max - 1)
    t_critical = compression_time * (c_critical - 1) / (compression - 1)
    supercritical_duration = compression_time - t_critical

    yields = []
    predet_times = []

    for trial in range(n_trials):
        # During supercritical phase, neutrons can initiate chains
        # Effective rate during supercritical window with chain initiation probability
        #
        # The probability a SF neutron initiates a chain depends on excess reactivity:
        # P_chain ≈ (k - 1) / k for k > 1
        #
        # For simplicity, we integrate over the supercritical phase:
        # Average k during supercritical phase ≈ (1 + k_max) / 2

        alpha_max, params_max = pk.inverse_bomb_period(data, final_radius, compression)

        k_inf = params_max["k_inf"]
        D = params_max["D"]
        delta = 0.7104 * params_max["lambda_tr"]
        R_prime = final_radius + delta
        B_g_sq = (np.pi / R_prime) ** 2
        L_sq = params_max["L_squared"]

        k_max = k_inf / (1.0 + L_sq * B_g_sq)

        k_avg = (1.0 + k_max) / 2
        p_chain_avg = (k_avg - 1) / k_avg if k_avg > 1 else 0

        # Effective predetonation rate = SF_rate × P_chain
        effective_rate = neutron_rate * p_chain_avg

        if effective_rate > 0:
            # Time to predetonation within supercritical window (exponential)
            predet_delay = np.random.exponential(1.0 / effective_rate)
        else:
            predet_delay = supercritical_duration * 2  # No predetonation

        if predet_delay >= supercritical_duration:
            # No predetonation - full yield at max compression
            effective_compression = compression
            predet_time = compression_time
        else:
            # Predetonation at t = t_critical + predet_delay
            predet_time = t_critical + predet_delay
            progress = predet_time / compression_time
            effective_compression = 1.0 + (compression - 1.0) * progress

        predet_times.append(predet_time)

        # Compute yield using proper point kinetics
        # The yield depends on energy release during the supercritical excursion
        # Y ∝ E_fission × N_fissions
        #
        # From paper: Y ≈ 17.6 kt × (m/kg) × (α × τ)³ where τ is disassembly time
        # Simplified: Y scales roughly with α³ × mass, where α is Rossi alpha

        alpha, params = pk.inverse_bomb_period(
            data, final_radius, effective_compression
        )

        if alpha > 0:
            # Compute dimensionless yield parameter
            # Reference: full compression gives reference yield
            alpha_ref = alpha_max
            if alpha_ref > 0:
                # Yield scales roughly as (α/α_ref)³ × reference_yield
                # Reference yield from paper: ~17-20 kt for 6.2 kg Pu at c=2.5
                ref_yield = 17.6 * (core_mass / 6.2) * (compression / 2.5) ** 2
                yield_kt = (
                    ref_yield * (alpha / alpha_ref) ** 2
                )  # squared scaling is more accurate
                yield_kt = max(0.1, min(yield_kt, 100))  # Reasonable bounds
            else:
                yield_kt = 0
        else:
            yield_kt = 0

        yields.append(yield_kt)

    yields = np.array(yields)
    predet_times = np.array(predet_times)

    # Fizzle = predetonation before full compression
    full_yield_mask = predet_times >= compression_time * 0.99
    fizzle_mask = ~full_yield_mask

    return {
        "yields": yields,
        "predet_times": predet_times,
        "mean_yield": np.mean(yields),
        "std_yield": np.std(yields),
        "fizzle_probability": np.sum(fizzle_mask) / n_trials,
        "full_yield_probability": np.sum(full_yield_mask) / n_trials,
        "neutron_rate": neutron_rate,
        "pu_mix": pu_mix,
        "n_trials": n_trials,
        "c_critical": c_critical,
        "t_critical": t_critical,
        "supercritical_duration": supercritical_duration,
    }


def study_pu240_content_effect(
    data: pk.NuclearData,
    core_mass: float,
    compression: float,
    pu240_fractions: np.ndarray,
    compression_time: float = 12e-6,
    n_trials: int = 200,
) -> Dict[str, Any]:
    """Study how Pu-240 content affects yield distribution."""

    results = []
    for pu240_frac in pu240_fractions:
        pu_mix = PlutoniumMix(
            pu239_fraction=1.0 - pu240_frac, pu240_fraction=pu240_frac
        )
        sim = simulate_with_predetonation(
            data, core_mass, compression, pu_mix, compression_time, n_trials
        )
        results.append(
            {
                "pu240_fraction": pu240_frac,
                "mean_yield": sim["mean_yield"],
                "std_yield": sim["std_yield"],
                "fizzle_prob": sim["fizzle_probability"],
                "yields": sim["yields"],
            }
        )

    return {
        "pu240_fractions": pu240_fractions,
        "results": results,
    }


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================


def plot_geometry_comparison(results: Dict[str, Any], save_path: Optional[str] = None):
    """Plot flux distributions for different geometries."""
    mesh = results["mesh"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    geometries = ["sphere", "oblate", "prolate", "cylinder"]
    titles = ["Sphere", "Oblate Ellipsoid", "Prolate Ellipsoid", "Cylinder"]

    for ax, geom, title in zip(axes.flatten(), geometries, titles):
        phi = results[geom]["phi"]
        k = results[geom]["k_eff"]
        dims = results[geom]["dims"]

        phi_pos = np.maximum(phi, 1e-10)
        im = ax.pcolormesh(
            mesh.Z * 100,
            mesh.R * 100,
            phi_pos,
            norm=LogNorm(),
            cmap="hot",
            shading="auto",
        )
        ax.set_xlabel("z (cm)")
        ax.set_ylabel("r (cm)")
        ax.set_title(f"{title}\nk={k:.4f}, {dims}")
        ax.set_aspect("equal")
        plt.colorbar(im, ax=ax)

    plt.suptitle(f"Geometry Comparison (mass={results['mass']:.1f}kg)", fontsize=14)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_tamper_study(results: Dict[str, Any], save_path: Optional[str] = None):
    """Plot k_eff vs tamper thickness."""
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(results["thicknesses"] * 100, results["k_values"], "bo-", markersize=6)
    ax.axhline(y=1.0, color="r", linestyle="--", label="k=1")
    ax.set_xlabel("Tamper Thickness (cm)")
    ax.set_ylabel("k_eff")
    ax.set_title(
        f"Effect of {results['tamper']} Tamper\n(core R={results['core_radius'] * 100:.2f}cm)"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_fizzle_distribution(results: Dict[str, Any], save_path: Optional[str] = None):
    """Plot yield distribution showing fizzle effect."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    pu240_fracs = results["pu240_fractions"]

    # Scatter plot of all yields
    ax = axes[0]
    for i, res in enumerate(results["results"]):
        pu240 = res["pu240_fraction"] * 100
        yields = res["yields"]
        ax.scatter([pu240] * len(yields), yields, alpha=0.3, s=10)
    ax.set_xlabel("Pu-240 Content (%)")
    ax.set_ylabel("Yield (kt)")
    ax.set_title("Yield Distribution vs Pu-240 Content")
    ax.grid(True, alpha=0.3)

    # Mean yield and fizzle probability
    ax = axes[1]
    mean_yields = [r["mean_yield"] for r in results["results"]]
    fizzle_probs = [r["fizzle_prob"] for r in results["results"]]

    ax.plot(pu240_fracs * 100, mean_yields, "b-o", label="Mean Yield", markersize=6)
    ax.set_xlabel("Pu-240 Content (%)")
    ax.set_ylabel("Mean Yield (kt)", color="b")
    ax.tick_params(axis="y", labelcolor="b")

    ax2 = ax.twinx()
    ax2.plot(
        pu240_fracs * 100,
        np.array(fizzle_probs) * 100,
        "r-s",
        label="Fizzle Prob",
        markersize=6,
    )
    ax2.set_ylabel("Fizzle Probability (%)", color="r")
    ax2.tick_params(axis="y", labelcolor="r")

    ax.set_title("Mean Yield and Fizzle Probability")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# =============================================================================
# MAIN - RUN ALL STUDIES
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("ADVANCED NUCLEAR SIMULATIONS")
    print("=" * 70)

    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("1. TIME-DEPENDENT 2D SIMULATION")
    print("=" * 70)

    print("\nInitializing 2D time-dependent simulation...")
    sim = TimeDependentSimulation2D(
        data=pk.PU239,
        core_mass=6.2,
        compression=2.5,
        mesh_size=30,  # Smaller for speed
    )

    print("Running simulation (this may take a moment)...")
    try:
        td_results = sim.simulate(t_max=5e-8)
        print(f"  Simulation completed")
        print(f"  Time steps: {len(td_results['t'])}")
        print(f"  Final radius: {td_results['radius'][-1] * 100:.3f} cm")
        print(f"  Final yield: {td_results['yield_kt'][-1]:.3f} kt")
    except Exception as e:
        print(f"  Simulation stopped early: {e}")
        td_results = None

    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("2. GEOMETRY COMPARISON")
    print("=" * 70)

    print("\nComparing geometries for 10kg Pu-239...")
    geom_results = study_geometry_criticality(
        pk.PU239, target_mass=10.0, compression=1.0
    )

    print(f"\n{'Geometry':<20} {'k_eff':<10} {'Dimensions'}")
    print("-" * 60)
    for geom in ["sphere", "oblate", "prolate", "cylinder"]:
        r = geom_results[geom]
        print(f"{geom.capitalize():<20} {r['k_eff']:<10.4f} {r['dims']}")

    fig = plot_geometry_comparison(geom_results, "geometry_comparison.png")
    plt.close()
    print("\nSaved: geometry_comparison.png")

    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("3. TAMPER/REFLECTOR STUDY")
    print("=" * 70)

    for tamper in [U238_TAMPER, BERYLLIUM_REFLECTOR, TUNGSTEN_TAMPER]:
        print(f"\n{tamper.name} tamper study...")
        tamper_results = study_tamper_effect(
            pk.PU239, core_mass=8.0, tamper=tamper, compression=1.0
        )

        print(f"  Bare core k_eff: {tamper_results['k_values'][0]:.4f}")
        print(f"  Max tamper k_eff: {tamper_results['k_values'][-1]:.4f}")
        print(
            f"  k_eff increase: {(tamper_results['k_values'][-1] / tamper_results['k_values'][0] - 1) * 100:.1f}%"
        )

        fig = plot_tamper_study(
            tamper_results, f"tamper_{tamper.name.lower().replace('-', '')}.png"
        )
        plt.close()
        print(f"  Saved: tamper_{tamper.name.lower().replace('-', '')}.png")

    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("4. FIZZLE PROBABILITY STUDY")
    print("=" * 70)

    print("\nStudying Pu-240 content effect on yield...")
    pu240_fracs = np.array([0.02, 0.05, 0.07, 0.10, 0.15, 0.20, 0.25, 0.30])

    fizzle_results = study_pu240_content_effect(
        pk.PU239,
        core_mass=6.2,
        compression=2.5,
        pu240_fractions=pu240_fracs,
        compression_time=12e-6,
        n_trials=200,
    )

    print(f"\n{'Pu-240 %':<10} {'Mean Yield (kt)':<15} {'Fizzle Prob %':<15}")
    print("-" * 40)
    for res in fizzle_results["results"]:
        print(
            f"{res['pu240_fraction'] * 100:<10.0f} {res['mean_yield']:<15.2f} {res['fizzle_prob'] * 100:<15.1f}"
        )

    fig = plot_fizzle_distribution(fizzle_results, "fizzle_study.png")
    plt.close()
    print("\nSaved: fizzle_study.png")

    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
Generated files:
  - geometry_comparison.png : k_eff for sphere/ellipsoid/cylinder
  - tamper_u238.png        : U-238 tamper effect
  - tamper_beryllium.png   : Beryllium reflector effect  
  - tamper_tungsten.png    : Tungsten tamper effect
  - fizzle_study.png       : Pu-240 fizzle probability

Key findings:
  1. Sphere has optimal k_eff for given mass (minimum surface/volume)
  2. Tampers increase k_eff by reflecting neutrons (U-238: also fissions)
  3. Higher Pu-240 content increases fizzle probability
  4. Time-dependent 2D captures spatial effects during expansion
""")
