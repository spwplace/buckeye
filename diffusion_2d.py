"""
2D Axisymmetric Neutron Diffusion - Finite Difference Implementation
Solves: ∂n/∂t = D∇²Φ + (νΣf - Σa)Φ
In (r,z): ∇² = ∂²/∂r² + (1/r)∂/∂r + ∂²/∂z²
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from dataclasses import dataclass
from typing import Dict, Any, Tuple, Optional
import point_kinetics as pk


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
        return i == 0 or i == self.nr - 1 or j == 0 or j == self.nz - 1


@dataclass
class Material2D:
    D: np.ndarray
    Sigma_f: np.ndarray
    Sigma_a: np.ndarray
    nu: float


def create_sphere_material(
    mesh: Mesh2D, data: pk.NuclearData, radius: float, compression: float
) -> Material2D:
    params = pk.compute_diffusion_params(data, compression)

    D = np.zeros((mesh.nr, mesh.nz))
    Sigma_f = np.zeros((mesh.nr, mesh.nz))
    Sigma_a = np.zeros((mesh.nr, mesh.nz))

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            r_dist = np.sqrt(mesh.R[i, j] ** 2 + mesh.Z[i, j] ** 2)
            if r_dist <= radius:
                D[i, j] = params["D"]
                Sigma_f[i, j] = params["Sigma_f"]
                Sigma_a[i, j] = params["Sigma_a"]
            else:
                D[i, j] = 1e-6
                Sigma_f[i, j] = 0
                Sigma_a[i, j] = 1e-6

    return Material2D(D=D, Sigma_f=Sigma_f, Sigma_a=Sigma_a, nu=data.nu_prompt)


def compute_k_effective(
    mesh: Mesh2D, mat: Material2D, tol: float = 1e-6, max_iter: int = 500
) -> Tuple[float, np.ndarray]:
    """
    k_eff via power iteration: LΦ = (1/k)FΦ
    L = -D∇² + Σa (leakage + absorption)
    F = νΣf (fission production)
    """
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

            # ∂²Φ/∂r² + (1/r)∂Φ/∂r  (axisymmetric Laplacian)
            if r > 1e-10:
                coeff_r_plus = D_local / dr**2 + D_local / (2 * r * dr)
                coeff_r_minus = D_local / dr**2 - D_local / (2 * r * dr)
                coeff_r_center = -2 * D_local / dr**2
            else:
                # r=0 singularity: L'Hopital gives (1/r)∂Φ/∂r → ∂²Φ/∂r²
                coeff_r_plus = 2 * D_local / dr**2
                coeff_r_minus = 2 * D_local / dr**2
                coeff_r_center = -4 * D_local / dr**2

            coeff_z = D_local / dz**2

            # L = -D∇² + Σa
            L[idx, mesh.node_index(i + 1, j)] = -coeff_r_plus
            L[idx, mesh.node_index(i - 1, j)] = -coeff_r_minus
            L[idx, mesh.node_index(i, j + 1)] = -coeff_z
            L[idx, mesh.node_index(i, j - 1)] = -coeff_z
            L[idx, idx] = -coeff_r_center + 2 * coeff_z + Sigma_a_local

            Fmat[idx, idx] = mat.nu * Sigma_f_local

    L_csr = csr_matrix(L)
    F_csr = csr_matrix(Fmat)

    phi = np.ones(n) * 0.1
    for i in range(mesh.nr):
        for j in range(mesh.nz):
            if mesh.is_boundary(i, j):
                phi[mesh.node_index(i, j)] = 0.0

    k = 1.0

    for iteration in range(max_iter):
        source = F_csr @ phi
        fission_rate = np.sum(source)

        if fission_rate < 1e-30:
            return 0.0, phi.reshape(mesh.nr, mesh.nz)

        phi_new = spsolve(L_csr, source / k)
        phi_new = np.maximum(phi_new, 0)

        new_fission_rate = np.sum(F_csr @ phi_new)
        if new_fission_rate > 1e-30:
            k_new = k * new_fission_rate / fission_rate
        else:
            k_new = 0.0

        norm = np.linalg.norm(phi_new)
        if norm > 1e-30:
            phi_new = phi_new / norm

        if abs(k_new - k) / max(k, 1e-10) < tol:
            return k_new, phi_new.reshape(mesh.nr, mesh.nz)

        k = k_new
        phi = phi_new

    return k, phi.reshape(mesh.nr, mesh.nz)


def compute_alpha_2d(mesh: Mesh2D, mat: Material2D) -> Tuple[float, np.ndarray]:
    """
    Dominant eigenvalue of dΦ/dt = (D∇² + νΣf - Σa)Φ
    Uses inverse power iteration on (A + F - αI).
    """
    n = mesh.n_nodes
    A = lil_matrix((n, n))

    dr, dz = mesh.dr, mesh.dz
    v = pk.NEUTRON_VELOCITY

    for i in range(mesh.nr):
        for j in range(mesh.nz):
            idx = mesh.node_index(i, j)
            r = mesh.r[i]

            if mesh.is_boundary(i, j):
                A[idx, idx] = -1e10
                continue

            D_local = mat.D[i, j]
            Sigma_a_local = mat.Sigma_a[i, j]
            Sigma_f_local = mat.Sigma_f[i, j]

            if r > 1e-10:
                coeff_r_plus = D_local / dr**2 + D_local / (2 * r * dr)
                coeff_r_minus = D_local / dr**2 - D_local / (2 * r * dr)
                coeff_r_center = -2 * D_local / dr**2
            else:
                coeff_r_plus = 2 * D_local / dr**2
                coeff_r_minus = 2 * D_local / dr**2
                coeff_r_center = -4 * D_local / dr**2

            coeff_z = D_local / dz**2

            # A = v*(D∇² + νΣf - Σa)
            A[idx, mesh.node_index(i + 1, j)] = v * coeff_r_plus
            A[idx, mesh.node_index(i - 1, j)] = v * coeff_r_minus
            A[idx, mesh.node_index(i, j + 1)] = v * coeff_z
            A[idx, mesh.node_index(i, j - 1)] = v * coeff_z
            A[idx, idx] = v * (
                coeff_r_center + 2 * coeff_z + mat.nu * Sigma_f_local - Sigma_a_local
            )

    A_csr = csr_matrix(A)

    phi = np.ones(n)
    for i in range(mesh.nr):
        for j in range(mesh.nz):
            if mesh.is_boundary(i, j):
                phi[mesh.node_index(i, j)] = 0.0

    alpha_est = 0.0

    for _ in range(200):
        phi_new = A_csr @ phi

        inner = np.inner(phi, phi_new)
        norm_phi = np.inner(phi, phi)
        if norm_phi > 1e-30:
            alpha_new = inner / norm_phi
        else:
            alpha_new = 0.0

        norm = np.linalg.norm(phi_new)
        if norm > 1e-30:
            phi = phi_new / norm

        if abs(alpha_new - alpha_est) < 1e-3 * abs(alpha_est) + 1e3:
            alpha_est = alpha_new
            break
        alpha_est = alpha_new

    return alpha_est, phi.reshape(mesh.nr, mesh.nz)


def plot_flux_distribution(mesh: Mesh2D, phi: np.ndarray, title: str = "Neutron Flux"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    phi_positive = np.maximum(phi, 1e-10)

    im = axes[0].pcolormesh(
        mesh.Z * 100,
        mesh.R * 100,
        phi_positive,
        norm=LogNorm(),
        cmap="hot",
        shading="auto",
    )
    axes[0].set_xlabel("z (cm)")
    axes[0].set_ylabel("r (cm)")
    axes[0].set_title(f"{title}")
    axes[0].set_aspect("equal")
    plt.colorbar(im, ax=axes[0])

    mid_z = mesh.nz // 2
    axes[1].plot(mesh.r * 100, phi[:, mid_z], "b-", linewidth=2)
    axes[1].set_xlabel("r (cm)")
    axes[1].set_ylabel("Flux")
    axes[1].set_title("Radial Profile at z=0")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    print("=" * 60)
    print("2D Axisymmetric Neutron Diffusion")
    print("=" * 60)

    mesh = Mesh2D(nr=60, nz=60, R_max=0.12, Z_max=0.24)

    print(f"\nMesh: {mesh.nr}x{mesh.nz} nodes")
    print(
        f"Domain: R=[0, {mesh.R_max * 100:.1f}cm], Z=[{-mesh.Z_max * 100 / 2:.1f}, {mesh.Z_max * 100 / 2:.1f}cm]"
    )

    print("\n" + "-" * 40)
    print("Critical radius scan (Pu-239, uncompressed)")
    print("-" * 40)

    mass_c, R_c, params = pk.compute_critical_mass(pk.PU239, compression=1.0)
    if R_c is not None:
        print(f"Point kinetics critical radius: {R_c * 100:.3f} cm")

        test_radii = np.linspace(R_c * 0.7, R_c * 1.3, 7)

        print(
            f"\n{'Radius (cm)':<12} {'k_eff (2D)':<12} {'alpha_point':<15} {'alpha_2d':<15}"
        )
        print("-" * 54)

        for radius in test_radii:
            mat = create_sphere_material(mesh, pk.PU239, radius, compression=1.0)
            k_eff, phi = compute_k_effective(mesh, mat)
            alpha_point, _ = pk.inverse_bomb_period(pk.PU239, radius, compression=1.0)
            alpha_2d, phi_alpha = compute_alpha_2d(mesh, mat)
            print(
                f"{radius * 100:<12.2f} {k_eff:<12.4f} {alpha_point:<15.2e} {alpha_2d:<15.2e}"
            )

    print("\n" + "-" * 40)
    print("Compressed core (6.2kg Pu-239)")
    print("-" * 40)

    core_mass = 6.2

    print(
        f"\n{'Compression':<12} {'Radius (cm)':<12} {'k_eff (2D)':<12} {'alpha_point':<15}"
    )
    print("-" * 51)

    for c in [1.5, 2.0, 2.5, 3.0]:
        compressed_density = pk.PU239.density_kg_m3 * c
        volume = core_mass / compressed_density
        radius = (3 * volume / (4 * np.pi)) ** (1 / 3)

        mat = create_sphere_material(mesh, pk.PU239, radius, compression=c)
        k_eff, phi = compute_k_effective(mesh, mat)
        alpha_point, _ = pk.inverse_bomb_period(pk.PU239, radius, compression=c)

        print(f"{c:<12.1f} {radius * 100:<12.2f} {k_eff:<12.4f} {alpha_point:<15.2e}")

    print("\n" + "-" * 40)
    print("Flux distribution plot")
    print("-" * 40)

    c = 2.5
    compressed_density = pk.PU239.density_kg_m3 * c
    volume = core_mass / compressed_density
    radius = (3 * volume / (4 * np.pi)) ** (1 / 3)

    mat = create_sphere_material(mesh, pk.PU239, radius, compression=c)
    k_eff, phi = compute_k_effective(mesh, mat)

    fig = plot_flux_distribution(mesh, phi, f"Pu-239 Eigenflux (c={c}, k={k_eff:.3f})")
    plt.savefig("flux_2d.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved: flux_2d.png")

    print("\n" + "-" * 40)
    print("k_eff vs radius")
    print("-" * 40)

    radii = np.linspace(0.025, 0.08, 15)
    k_values = []
    alpha_values = []

    for r in radii:
        mat = create_sphere_material(mesh, pk.PU239, r, compression=1.0)
        k, _ = compute_k_effective(mesh, mat)
        alpha, _ = pk.inverse_bomb_period(pk.PU239, r, compression=1.0)
        k_values.append(k)
        alpha_values.append(alpha)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(radii * 100, k_values, "b-o", markersize=4)
    ax1.axhline(y=1.0, color="r", linestyle="--", label="k=1")
    ax1.set_xlabel("Radius (cm)")
    ax1.set_ylabel("k_eff")
    ax1.set_title("2D FEM k_eff vs Radius")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(radii * 100, np.array(alpha_values) / 1e8, "g-o", markersize=4)
    ax2.axhline(y=0, color="r", linestyle="--", label="alpha=0")
    ax2.set_xlabel("Radius (cm)")
    ax2.set_ylabel("Alpha (10^8 /s)")
    ax2.set_title("Point Kinetics Alpha vs Radius")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("criticality_transition.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved: criticality_transition.png")
