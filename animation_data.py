"""
Animation Data Generator - pre-computes simulation data for manim animations.
Units: SI (meters, seconds, Kelvin, Pascals, Joules) unless noted.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional

import nuclear_physics as nuc
import spatial_neutronics as sn


@dataclass
class ImplosionData:
    t: np.ndarray  # seconds
    R_core: np.ndarray  # Pu-239 core radius (m)
    R_tamper: np.ndarray  # U-238 tamper radius (m)
    R_pusher: np.ndarray  # Al pusher radius (m)
    R_outer: np.ndarray  # HE/casing boundary (m)
    compression: np.ndarray  # ρ/ρ₀
    v_core: np.ndarray  # core surface velocity (m/s)
    i_critical: int  # index when k=1 crossed
    i_peak_compression: int  # index of max compression
    core_mass: float
    tamper_mass: float
    pusher_mass: float
    outer_mass: float


@dataclass
class ExcursionData:
    t: np.ndarray  # seconds
    R: np.ndarray  # radius (m)
    v: np.ndarray  # velocity (m/s)
    N: np.ndarray  # neutron count
    alpha: np.ndarray  # Rossi alpha (1/s)
    T: np.ndarray  # temperature (K)
    P: np.ndarray  # pressure (Pa)
    E_fission: np.ndarray  # fission energy (J)
    E_kinetic: np.ndarray  # kinetic energy (J)
    E_internal: np.ndarray  # internal energy (J)
    yield_kt: np.ndarray  # yield (kilotons TNT)
    compression: np.ndarray  # ρ/ρ₀
    i_peak_alpha: int  # index of max reactivity
    i_peak_power: int  # index of max fission rate
    i_subcritical: int  # index when α crosses zero


@dataclass
class FluxFieldData:
    Nr: int
    Nz: int
    r: np.ndarray  # radial coords (m)
    z: np.ndarray  # axial coords (m)
    R_max: float
    Z_max: float
    phi: np.ndarray  # flux field (Nr × Nz)
    material_mask: np.ndarray  # boolean mask for fissile region
    k_eff: float
    sphere_radius: float
    compression: float


@dataclass
class FVMCellData:
    r_centers: np.ndarray
    z_centers: np.ndarray
    r_faces: np.ndarray
    z_faces: np.ndarray
    volumes: np.ndarray
    D: np.ndarray  # diffusion coefficient
    Sigma_f: np.ndarray  # fission cross-section
    Sigma_a: np.ndarray  # absorption cross-section
    J_r: np.ndarray  # radial current at r-faces
    J_z: np.ndarray  # axial current at z-faces
    source: np.ndarray  # fission source (νΣ_f × φ)


@dataclass
class TimeEvolutionFlux:
    t: np.ndarray
    phi_sequence: np.ndarray  # shape: (n_times, Nr, Nz)
    R_sequence: np.ndarray  # sphere radius at each time
    r: np.ndarray
    z: np.ndarray


def generate_implosion_data(
    core_mass: float = 6.2,
    tamper_mass: float = 108.0,
    pusher_mass: float = 130.0,
    outer_mass: float = 4430.0,
    peak_compression: float = 2.5,
    compression_time: float = 12e-6,
    n_points: int = 500,
) -> ImplosionData:
    shell = nuc.MassShell(core_mass, tamper_mass, pusher_mass, outer_mass)

    rho_Pu = nuc.PU239.rho_0
    V_core = core_mass / rho_Pu
    R_core_0 = (3 * V_core / (4 * np.pi)) ** (1 / 3)

    t = np.linspace(0, compression_time, n_points)

    t_norm = t / compression_time
    compression = 1.0 + (peak_compression - 1.0) * (
        0.5 * (1 + np.tanh(6 * (t_norm - 0.5)))
    )

    scale = compression ** (-1 / 3)

    R_core = R_core_0 * scale
    R_tamper = shell.R_tamper_ref * scale
    R_pusher = shell.R_pusher_ref * scale
    R_outer = shell.R_outer_ref * np.ones_like(t)

    v_core = np.gradient(R_core, t)

    i_critical = 0
    for i, c in enumerate(compression):
        try:
            m_c, R_c, params = nuc.compute_critical_mass(nuc.PU239, c)
            if core_mass >= m_c:
                i_critical = i
                break
        except ValueError:
            continue

    i_peak = int(np.argmax(compression))

    return ImplosionData(
        t=t,
        R_core=R_core,
        R_tamper=R_tamper,
        R_pusher=R_pusher,
        R_outer=R_outer,
        compression=compression,
        v_core=v_core,
        i_critical=i_critical,
        i_peak_compression=i_peak,
        core_mass=core_mass,
        tamper_mass=tamper_mass,
        pusher_mass=pusher_mass,
        outer_mass=outer_mass,
    )


def generate_excursion_data(
    core_mass: float = 6.2,
    compression: float = 2.5,
    N_0: float = 1e8,
    t_max: float = 2e-6,
) -> ExcursionData:
    sim = nuc.SupercriticalSimulation(
        data=nuc.PU239,
        core_mass=core_mass,
        compression=compression,
        use_full_mass_shell=False,
    )

    results = sim.simulate(N_0=N_0, t_max=t_max)

    if "error" in results:
        raise ValueError(f"Simulation failed: {results['error']}")

    i_peak_alpha = int(np.argmax(results["alpha"]))

    power = np.gradient(results["E_fission"], results["t"])
    i_peak_power = int(np.argmax(power))

    alpha = results["alpha"]
    i_subcritical = len(alpha) - 1
    for i in range(len(alpha) - 1):
        if alpha[i] > 0 and alpha[i + 1] <= 0:
            i_subcritical = i
            break

    return ExcursionData(
        t=results["t"],
        R=results["R"],
        v=results["v"],
        N=results["N"],
        alpha=results["alpha"],
        T=results["T"],
        P=results["P"],
        E_fission=results["E_fission"],
        E_kinetic=results["E_kinetic"],
        E_internal=results["E_internal"],
        yield_kt=results["yield_kt"],
        compression=results["compression"],
        i_peak_alpha=i_peak_alpha,
        i_peak_power=i_peak_power,
        i_subcritical=i_subcritical,
    )


def generate_flux_field(
    compression: float = 1.0, Nr: int = 60, Nz: int = 120, domain_factor: float = 3.0
) -> FluxFieldData:
    m_c, R_c, params = nuc.compute_critical_mass(nuc.PU239, compression)

    core_mass = 6.2
    rho = nuc.PU239.rho_0 * compression
    V = core_mass / rho
    sphere_radius = (3 * V / (4 * np.pi)) ** (1 / 3)

    domain_size = sphere_radius * domain_factor
    mesh = sn.CylindricalMesh(Nr=Nr, Nz=Nz, R_max=domain_size, Z_max=2 * domain_size)

    mat = sn.create_homogeneous_sphere(mesh, nuc.PU239, sphere_radius, compression)
    k_eff, phi = sn.compute_k_effective(mesh, mat)

    return FluxFieldData(
        Nr=Nr,
        Nz=Nz,
        r=mesh.r,
        z=mesh.z,
        R_max=mesh.R_max,
        Z_max=mesh.Z_max,
        phi=phi,
        material_mask=mat.volume_mask,
        k_eff=k_eff,
        sphere_radius=sphere_radius,
        compression=compression,
    )


def generate_fvm_cell_data(
    compression: float = 2.5, Nr: int = 30, Nz: int = 60
) -> FVMCellData:
    core_mass = 6.2
    rho = nuc.PU239.rho_0 * compression
    V = core_mass / rho
    sphere_radius = (3 * V / (4 * np.pi)) ** (1 / 3)

    domain_size = sphere_radius * 2.5
    mesh = sn.CylindricalMesh(Nr=Nr, Nz=Nz, R_max=domain_size, Z_max=2 * domain_size)

    mat = sn.create_homogeneous_sphere(mesh, nuc.PU239, sphere_radius, compression)
    k_eff, phi = sn.compute_k_effective(mesh, mat)

    dr, dz = mesh.dr, mesh.dz

    J_r = np.zeros((Nr + 1, Nz))
    for i in range(1, Nr):
        for j in range(Nz):
            D_avg = 0.5 * (mat.D[i - 1, j] + mat.D[i, j])
            J_r[i, j] = -D_avg * (phi[i, j] - phi[i - 1, j]) / dr

    J_z = np.zeros((Nr, Nz + 1))
    for i in range(Nr):
        for j in range(1, Nz):
            D_avg = 0.5 * (mat.D[i, j - 1] + mat.D[i, j])
            J_z[i, j] = -D_avg * (phi[i, j] - phi[i, j - 1]) / dz

    source = mat.nu_Sigma_f * phi

    return FVMCellData(
        r_centers=mesh.r,
        z_centers=mesh.z,
        r_faces=mesh.r_face,
        z_faces=mesh.z_face,
        volumes=mesh.V,
        D=mat.D,
        Sigma_f=mat.Sigma_f,
        Sigma_a=mat.Sigma_a,
        J_r=J_r,
        J_z=J_z,
        source=source,
    )


def generate_time_evolution_flux(
    compression_start: float = 2.0,
    compression_end: float = 1.5,
    n_frames: int = 50,
    Nr: int = 40,
    Nz: int = 80,
) -> TimeEvolutionFlux:
    compressions = np.linspace(compression_start, compression_end, n_frames)
    t = np.linspace(0, 1e-6, n_frames)

    phi_sequence = []
    R_sequence = []

    core_mass = 6.2

    for c in compressions:
        rho = nuc.PU239.rho_0 * c
        V = core_mass / rho
        R = (3 * V / (4 * np.pi)) ** (1 / 3)
        R_sequence.append(R)

        domain_size = R * 2.5
        mesh = sn.CylindricalMesh(
            Nr=Nr, Nz=Nz, R_max=domain_size, Z_max=2 * domain_size
        )
        mat = sn.create_homogeneous_sphere(mesh, nuc.PU239, R, c)

        try:
            k_eff, phi = sn.compute_k_effective(mesh, mat)
            phi_max = np.max(phi)
            if phi_max > 0:
                phi = phi / phi_max
            if k_eff > 1:
                phi = phi * np.exp(2 * (k_eff - 1))
        except Exception:
            phi = np.zeros((Nr, Nz))

        phi_sequence.append(phi)

    mid_c = (compression_start + compression_end) / 2
    rho = nuc.PU239.rho_0 * mid_c
    V = core_mass / rho
    R = (3 * V / (4 * np.pi)) ** (1 / 3)
    domain_size = R * 2.5
    mesh = sn.CylindricalMesh(Nr=Nr, Nz=Nz, R_max=domain_size, Z_max=2 * domain_size)

    return TimeEvolutionFlux(
        t=t,
        phi_sequence=np.array(phi_sequence),
        R_sequence=np.array(R_sequence),
        r=mesh.r,
        z=mesh.z,
    )


@dataclass
class AnimationDataBundle:
    implosion: ImplosionData
    excursion: ExcursionData
    flux_field: FluxFieldData
    fvm_cells: FVMCellData
    time_evolution: TimeEvolutionFlux
    core_mass: float = 6.2
    material: str = "Pu-239"
    peak_compression: float = 2.5


def generate_all_animation_data(verbose: bool = True) -> AnimationDataBundle:
    if verbose:
        print("=" * 60)
        print("GENERATING ANIMATION DATA")
        print("=" * 60)

    if verbose:
        print("\n1. Generating implosion trajectory...")
    implosion = generate_implosion_data()
    if verbose:
        print(
            f"   {len(implosion.t)} time points, peak compression {implosion.compression.max():.2f}x"
        )

    if verbose:
        print("\n2. Running supercritical excursion simulation...")
    excursion = generate_excursion_data()
    if verbose:
        print(f"   Final yield: {excursion.yield_kt[-1]:.2f} kt")
        print(f"   Peak temperature: {excursion.T.max():.2e} K")

    if verbose:
        print("\n3. Computing 2D flux field...")
    flux_field = generate_flux_field(compression=2.5)
    if verbose:
        print(f"   k_eff = {flux_field.k_eff:.4f}")

    if verbose:
        print("\n4. Generating FVM cell data...")
    fvm_cells = generate_fvm_cell_data()
    if verbose:
        print(
            f"   {fvm_cells.r_centers.shape[0]} x {fvm_cells.z_centers.shape[0]} cells"
        )

    if verbose:
        print("\n5. Computing time-evolution flux sequence...")
    time_evolution = generate_time_evolution_flux(n_frames=30)
    if verbose:
        print(f"   {len(time_evolution.t)} frames")

    if verbose:
        print("\n" + "=" * 60)
        print("DATA GENERATION COMPLETE")
        print("=" * 60)

    return AnimationDataBundle(
        implosion=implosion,
        excursion=excursion,
        flux_field=flux_field,
        fvm_cells=fvm_cells,
        time_evolution=time_evolution,
    )


_cached_data: Optional[AnimationDataBundle] = None


def get_animation_data(regenerate: bool = False) -> AnimationDataBundle:
    global _cached_data
    if _cached_data is None or regenerate:
        _cached_data = generate_all_animation_data()
    return _cached_data


if __name__ == "__main__":
    data = generate_all_animation_data(verbose=True)

    print("\n--- Implosion Data ---")
    print(
        f"Critical crossing at t = {data.implosion.t[data.implosion.i_critical] * 1e6:.2f} μs"
    )
    print(
        f"Core radius: {data.implosion.R_core[0] * 100:.2f} → {data.implosion.R_core[-1] * 100:.2f} cm"
    )

    print("\n--- Excursion Data ---")
    print(f"Duration: {data.excursion.t[-1] * 1e6:.3f} μs")
    print(f"Peak α: {data.excursion.alpha[data.excursion.i_peak_alpha]:.2e} /s")

    print("\n--- Flux Field ---")
    print(f"Grid: {data.flux_field.Nr} x {data.flux_field.Nz}")
    print(f"k_eff: {data.flux_field.k_eff:.4f}")
