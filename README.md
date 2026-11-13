# CFD Solver Apps

This repository contains two MATLAB apps:

- **unsteady_diffusion_app** - a graphical interface for the 2D unsteady diffusion solver.
- **cfd_app** - a combined app for the 2D unsteady diffusion solver and the 2D Navier-Stokes lid-driven cavity solver.

## Recommended method: install the packaged apps

The easiest way to run the apps is to download the packaged MATLAB app files from the main page of this repository.

### unsteady_diffusion_app

1. Download `Unsteady_diffusion_app.mltbx` from the repository.
2. Open MATLAB.
3. Double-click the downloaded `.mltbx` file, or use **Apps > Install App** and select it.
4. Complete the installation.
5. Open the app from the MATLAB **Apps** gallery.

### cfd_app

1. Download `cfd_app.mltbx` from the repository.
2. Open MATLAB.
3. Double-click the downloaded `.mltbx` file, or use **Apps > Install App** and select it.
4. Complete the installation.
5. Open the app from the MATLAB **Apps** gallery.


## If the packaged app does not work

If installation of an `.mltbx` file is unsuccessful, the apps can also be run directly from their MATLAB source files in this ZIP.


## Running unsteady_diffusion_app from source

1. Open MATLAB.
2. Navigate to the extracted `unsteady_diffusion_app` folder and set it as the MATLAB **Current Folder**.
3. Make sure both files are present:
   - `unsteady_diffusion_app.m`
   - `solve_unsteady_diff.m`
4. In the MATLAB Command Window, run:

```matlab
unsteady_diffusion_app
```

The app window should open.

## Running cfd_app from source

1. Open MATLAB.
2. Navigate to the extracted `cfd_app` folder and set it as the MATLAB **Current Folder**.
3. Make sure the following files are present:
   - `cfd_app.m`
   - `solve_unsteady_diff.m`
   - `navier_stokes_2d.m`
4. In the MATLAB Command Window, run:

```matlab
cfd_app
```

The CFD app window should open.

## MATLAB path alternative

If MATLAB is not currently inside the app folder, add the required folder to the MATLAB path before running the app by right clicking on the folder and selecting "add to path".
