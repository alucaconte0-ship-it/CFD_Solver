# CFD Solver

MATLAB project developed for ME2610 to explore numerical solutions of partial differential equations using finite-difference methods.

The current project contains a **2D unsteady diffusion solver** and a graphical MATLAB app for setting up and visualising diffusion simulations. A separate 2D incompressible Navier-Stokes lid-driven cavity solver is also included as an extension.

## 2D unsteady diffusion solver

`solve_unsteady_diff.m` solves the two-dimensional unsteady diffusion equation

```text
du/dt = nu (d2u/dx2 + d2u/dy2)
```

using an explicit Forward Euler time discretisation and second-order central finite differences in space. The spatial discretisation gives the standard five-point stencil.

The solver supports:

- rectangular domains defined by `Lx` and `Ly`
- independently specified grid sizes `Mx` and `My`
- constant Dirichlet boundary values on all four sides
- optional periodic boundary conditions in the x and y directions
- arbitrary `Mx x My` initial-condition arrays
- optional plotting while the simulation runs

The main function is:

```matlab
u = solve_unsteady_diff(uinitial,Lx,Ly,Mx,My,nu,maxIter,dt, ...
    xperiodic,yperiodic,tBC,bBC,lBC,rBC,showplot)
```

The returned matrix `u` contains the scalar field after `maxIter` time steps.

## Unsteady diffusion app

`unsteady_diffusion_app.m` provides a graphical interface for the diffusion solver. It allows the user to change the domain dimensions, grid resolution, diffusivity, time step, number of iterations and boundary conditions without editing MATLAB code directly.

Three initial conditions are currently available:

- **Gaussian** - a smooth hotspot centred in the domain
- **Single hot cell** - one central grid point is assigned the initial amplitude
- **Uniform** - the whole domain starts at the same value

The app calls `solve_unsteady_diff.m` with plotting disabled inside the solver, then displays the final solution as a filled contour plot using MATLAB's `turbo` colormap.

## Repository structure

```text
CFD_Solver/
├── README.md
├── diffusion/
│   ├── solve_unsteady_diff.m
│   └── unsteady_diffusion_app.m
└── navier_stokes/
    └── navier_stokes_2d.m
```

## Running the diffusion app

Open MATLAB and make the `diffusion` folder the current folder, then run:

```matlab
unsteady_diffusion_app
```

Alternatively, add the project folders to the MATLAB path:

```matlab
addpath('diffusion')
addpath('navier_stokes')
```

The diffusion app only requires the files in the `diffusion` folder.

## Numerical note

Because the diffusion solver uses an explicit time-stepping scheme, the time step must be chosen small enough for numerical stability. For a uniform rectangular grid, the usual 2D stability condition is

```text
nu*dt*(1/dx^2 + 1/dy^2) <= 0.5
```

Further verification, validation and documentation will be added as the project develops.
