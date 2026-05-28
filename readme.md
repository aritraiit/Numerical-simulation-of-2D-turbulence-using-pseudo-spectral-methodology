# 2D Pseudo-Spectral Navier-Stokes Solver

This repository contains a high-order numerical framework for solving the 2D incompressible Navier-Stokes equations. The solver is built on the vorticity-streamfunction (ω-ψ) formulation and leverages a pseudo-spectral strategy. This approach combines the extreme spatial accuracy of Fourier spectral methods with the computational efficiency of evaluating nonlinear terms in physical space.

**Author:** Aritra Roy, former *S.N. Bhatt Memorial Excellence Fellow (2024)* at the International Centre for Theoretical Sciences (ICTS-TIFR), Bangalore.

---

## Technical Specifications & Mathematical Formulation

The solver models fluid flow in a doubly-periodic square domain (x, y) ∈ [0, 2π) × [0, 2π). 

The governing 2D advection-diffusion equation for vorticity is:

∂ω/∂t + u · ∇ω = ν ∇² ω

Where:
* ω = ∇ × u = -∇²ψ is the scalar vorticity field.
* ψ is the streamfunction, mapping the velocity components: u = ∂ψ/∂y and v = -∂ψ/∂x.
* ν is the kinematic viscosity parameter.

### Pseudo-Spectral Strategy
By transforming the equations into Fourier space, spatial derivatives and Laplacian inversions become highly efficient, local algebraic multiplications involving discrete wavenumbers k_x and k_y:
* Gradient: ∂/∂x → i k_x
* Laplacian: ∇² → -k² = -(k_x² + k_y²)

The nonlinear convective term (u · ∇ω) is evaluated point-by-point in **physical space** to avoid executing expensive O(N⁴) mathematical convolutions in frequency space.

### 2/3 Dealiasing Rule (Orszag Filter)
Because nonlinear terms are multiplied in the physical grid, high-frequency energy truncation errors (aliasing) are generated. To clear out these unresolvable grid modes and guarantee numerical stability, **Orszag’s 2/3 rule** is applied as a sharp filter in Fourier space:

ω̂(k_x, k_y) = 0 for √(k_x² + k_y²) > N/3

---

## Code Overview & Architecture

The script handles the simulation through several optimized computational steps:

1. **Grid Generation:** The physical mesh uses N_x × N_y points. Using `numpy.linspace` with `endpoint=False` omits the 2π boundary to satisfy periodic conditions without coordinate duplication.
2. **Fourier Wavenumbers:** Frequencies are mapped using `np.fft.fftfreq` and scaled to the physical domain. The Orszag filter is immediately applied.
3. **Precomputed Spectral Operators:** The Laplacian inverses and derivative multipliers are pre-allocated in memory to rapidly construct velocity fields from vorticity.
4. **Initialization:** The flow field is initialized directly in Fourier space. Injecting energy solely into low Fourier modes produces a perfectly smooth initial physical field without numerical shocks.
5. **Time Integration:** The PDE is flattened into a stiff ODE system (dΩ/dt = RHS(Ω)) and advanced through time using adaptive time-stepping (`scipy.integrate.solve_ivp` with RK45).

---

## Dependencies

* `numpy`
* `scipy`
* `matplotlib`
* `ipython` (for live notebook rendering)
* A local LaTeX installation (if `text.usetex` is set to `True` for rendering plots)

---

## Complete Source Code

The complete source code is included in this repository.

