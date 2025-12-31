# 2D Finite Volume Methods

Point kinetics treats the entire system as a single point. Real neutron distributions are spatially varying—higher in the center, zero at the boundary. The Python code includes 2D spatial neutronics using the **finite volume method** (FVM).

## Why 2D?

The nuclear football has axial symmetry—it's rotationally symmetric around its long axis. This means we can reduce the 3D problem to 2D in (r, z) cylindrical coordinates.

For a sphere, we could even reduce to 1D (radial only), but 2D allows us to:
- Study non-spherical geometries
- Resolve axial flux variations
- Validate against 3D Monte Carlo

## The Cylindrical Diffusion Equation

In cylindrical (r, z) coordinates, the diffusion equation becomes:

$$-\nabla \cdot (D\nabla\phi) + \Sigma_a\phi = \nu\Sigma_f\phi$$

where the Laplacian is:

$$\nabla^2\phi = \frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial\phi}{\partial r}\right) + \frac{\partial^2\phi}{\partial z^2}$$

The 1/r term creates a coordinate singularity at r = 0 (the axis), which requires special treatment.

## Finite Volume Discretization

We divide the domain into cells and integrate the diffusion equation over each cell. For cell (i, j) with volume V:

$$\int_V \left[-\nabla \cdot (D\nabla\phi) + \Sigma_a\phi - \nu\Sigma_f\phi\right] dV = 0$$

Using the divergence theorem:

$$-\oint_S D\nabla\phi \cdot d\mathbf{A} + V\Sigma_a\phi - V\nu\Sigma_f\phi = 0$$

The surface integral becomes flux through cell faces.

## Mesh Structure

The Python code uses a regular cylindrical mesh:

```python
@dataclass
class CylindricalMesh:
    Nr: int      # cells in r
    Nz: int      # cells in z
    R_max: float # outer radius
    Z_max: float # total height
    
    def __post_init__(self):
        self.dr = self.R_max / self.Nr
        self.dz = self.Z_max / self.Nz
        
        # Cell centers
        self.r = np.linspace(self.dr/2, self.R_max - self.dr/2, self.Nr)
        self.z = np.linspace(-self.Z_max/2 + self.dz/2, self.Z_max/2 - self.dz/2, self.Nz)
        
        # Cell volumes: V = 2πr dr dz
        self.V = 2 * np.pi * self.R * self.dr * self.dz
        self.V[0, :] = np.pi * (self.dr/2)**2 * self.dz  # axis cell
```

The axis cell (i = 0) has a different volume formula because it's a cylinder, not an annulus.

## Face Fluxes

The key to FVM is computing fluxes through cell faces. For the r-direction face between cells (i, j) and (i+1, j):

$$J_r = -D \frac{\partial\phi}{\partial r} \approx -D_{face} \frac{\phi_{i+1,j} - \phi_{i,j}}{\Delta r}$$

The face diffusion coefficient uses the **harmonic mean**:

$$D_{face} = \frac{2 D_i D_{i+1}}{D_i + D_{i+1}}$$

This ensures correct behavior at material interfaces.

## Boundary Conditions

### Axis of Symmetry (r = 0)
By symmetry: ∂φ/∂r = 0. No flux crosses the axis.

### Material-Vacuum Interface
Neutrons leak out at material boundaries. The extrapolated boundary condition:

$$\phi(r_{boundary} + \delta) = 0$$

where δ = 0.7104 λ_tr is the extrapolation distance.

This is implemented by adding a leakage term:

```python
def vacuum_leakage_coeff(D_mat, A, d_half):
    delta = 0.7104 * lambda_tr
    d_extrap = d_half + delta
    return D_mat * A / d_extrap
```

## Matrix Assembly

The discretized equations form a linear system:

$$\mathbf{L}\boldsymbol{\phi} = \frac{1}{k}\mathbf{F}\boldsymbol{\phi}$$

where:
- **L** is the loss matrix (leakage + absorption)
- **F** is the fission matrix

```python
def build_fvm_matrices(mesh, mat):
    L = lil_matrix((n, n))
    F = lil_matrix((n, n))
    
    for i in range(mesh.Nr):
        for j in range(mesh.Nz):
            idx = mesh.cell_index(i, j)
            V = mesh.V[i, j]
            
            # Leakage coefficients
            coeff_rp = D_face * A_r_plus / dr
            coeff_rm = D_face * A_r_minus / dr
            coeff_zp = D_face * A_z_plus / dz
            coeff_zm = D_face * A_z_minus / dz
            
            # Fill L matrix
            L[idx, idx] = (coeff_rp + coeff_rm + coeff_zp + coeff_zm)/V + Sigma_a
            L[idx, mesh.cell_index(i+1, j)] = -coeff_rp / V
            # ... other neighbors
            
            # Fill F matrix
            F[idx, idx] = nu_Sigma_f
    
    return csr_matrix(L), csr_matrix(F)
```

## Power Iteration

We solve the eigenvalue problem using power iteration:

```python
def power_iteration(L, F, tol=1e-7, max_iter=1000):
    phi = np.ones(n) / np.linalg.norm(np.ones(n))
    k = 1.0
    
    for iteration in range(max_iter):
        source = F @ phi
        total_source = np.sum(source)
        
        phi_new = spsolve(L, source / k)
        phi_new = np.maximum(phi_new, 0)
        
        new_total_source = np.sum(F @ phi_new)
        k_new = k * new_total_source / total_source
        
        phi_new = phi_new / np.linalg.norm(phi_new)
        
        if abs(k_new - k) / k < tol:
            return k_new, phi_new
        
        k = k_new
        phi = phi_new
    
    return k, phi
```

This converges to the dominant eigenvalue (k_eff) and eigenvector (flux shape).

## Mesh Convergence

The discretization introduces truncation error. We verify the code by checking mesh convergence:

```
N=20:  k_eff = 1.0234
N=30:  k_eff = 1.0156
N=40:  k_eff = 1.0118
N=50:  k_eff = 1.0098
N=60:  k_eff = 1.0085

Richardson extrapolated: k_eff = 1.0042
Convergence order: 1.9 (expected ~2 for 2nd order FVM)
```

The method converges at nearly the theoretical rate.

## Results for Pu-239 Critical Sphere

| Quantity | Point Kinetics | 2D FVM | Monte Carlo |
|----------|---------------|--------|-------------|
| k_eff at R_c | 1.000 | 1.004 | 1.057 |
| Critical radius | 4.9 cm | 5.0 cm | 6.4 cm |
| Critical mass | 10.0 kg | 10.4 kg | 16.3 kg |

The 2D FVM agrees well with point kinetics (both use diffusion theory) but differs from Monte Carlo due to the diffusion approximation.

## Flux Distribution

The calculated flux profile shows the expected shape:
- Maximum at center (r = 0, z = 0)
- Falls as sin(πr/R')/r in the radial direction
- Falls as sin(πz/Z') in the axial direction
- Zero at extrapolated boundary

## Limitations

The 2D FVM approach inherits diffusion theory limitations:
- Poor near vacuum boundaries
- Approximate for thin systems
- One-group energy treatment

But it adds:
- Spatial resolution
- Non-uniform material handling
- Flux shape calculation

## When to Use 2D FVM

**Good for:**
- Flux shape visualization
- Multi-region problems
- Material interface effects
- Verifying point kinetics

**Not for:**
- Accurate critical mass (use Monte Carlo)
- Complex geometries (use Monte Carlo)
- Near-boundary physics (use transport)

## Comparison Summary

| Method | Geometry | Energy | Speed | Accuracy |
|--------|----------|--------|-------|----------|
| Point kinetics | Sphere | 1 group | ~1 s | ~15% |
| 2D FVM | Any 2D | 1-2 groups | ~10 s | ~10% |
| Monte Carlo | Any 3D | Continuous | ~100 s | ~1-5% |

The 2D FVM provides intermediate fidelity between point kinetics and Monte Carlo, useful for design studies and physical understanding.
