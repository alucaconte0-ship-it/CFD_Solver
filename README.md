# CFD Solver App

This repository contains a MATLAB application for:

- 2D unsteady diffusion
- 2D incompressible Navier–Stokes lid-driven cavity flow

The app provides a graphical interface for entering simulation parameters, selecting the solver and visualising the numerical result.

## Recommended method: install the packaged app

The easiest way to run the program is to install the packaged MATLAB app.

1. Download `cfd_app.mltbx` from the repository.
2. Open MATLAB.
3. Double-click the downloaded `.mltbx` file, or use **Apps > Install App** and select it.
4. Complete the installation.
5. Open the app from the MATLAB **Apps** gallery.

## Running from source

If you are unable to install and run the packaged app, the application can also be run directly from the MATLAB source files.

1. Open MATLAB.
2. Navigate to the extracted `cfd_app` folder and set it as the MATLAB **Current Folder** or right click and select **Add to Path**.
3. Make sure the following files are present:
   - `cfd_app.m`
   - `solve_unsteady_diff.m`
   - `navier_stokes_2d.m`
4. In the MATLAB Command Window, run:

```matlab
cfd_app
```

The application window should open.

## Solver options

### 2D Diffusion

The diffusion mode solves the 2D unsteady diffusion equation using the `solve_unsteady_diff` function.

The user can define:

- domain dimensions
- grid resolution
- diffusivity
- timestep and number of timesteps
- fixed or periodic boundary conditions
- initial condition

### 2D Navier–Stokes

The Navier–Stokes mode solves a 2D incompressible lid-driven cavity flow using the `navier_stokes_2d` function.

The user can define:

- domain dimensions
- grid resolution
- kinematic viscosity
- timestep and number of timesteps
- fluid density
- lid velocity
- number of pressure iterations

The result is displayed as velocity magnitude with velocity vectors overlaid.
