# reflectomITER

`reflectomITER` is an AMReX-based finite-difference time-domain electromagnetic solver for Maxwell/plasma-wave simulations.  It supports MPI parallelism through AMReX, optional CUDA acceleration, optional IMAS/IMASH equilibrium loading, perfectly matched layers (PML), divergence cleaning, localized electromagnetic sources, and an auxiliary-differential-equation (ADE) plasma response model.

The code can be run in two main modes:

1. **Manual-domain mode**, where the user directly specifies the computational box and grid.
2. **IMAS/IMASH mode**, where the domain and plasma/background-field data are derived from an ITER/IMAS equilibrium.

---

## Contents

- [Main capabilities](#main-capabilities)
- [Code organization](#code-organization)
- [Theory and numerical model](#theory-and-numerical-model)
- [Dependencies](#dependencies)
- [Installing AMReX](#installing-amrex)
- [Building reflectomITER](#building-fdtd_amrex)
- [Running the code](#running-the-code)
- [Input-file structure](#input-file-structure)
- [Example 1: manual-domain run](#example-1-manual-domain-run)
- [Example 2: IMAS/IMASH run](#example-2-imasimash-run)
- [Output](#output)
- [Generating Doxygen documentation](#generating-doxygen-documentation)
- [Troubleshooting](#troubleshooting)

---

## Main capabilities

`reflectomITER` currently provides:

- electromagnetic FDTD updates for electric and magnetic fields,
- optional current-density coupling through an `EHJSolver`,
- localized sources deposited into `J`, `E`, or `H`,
- Gaussian, super-Gaussian, and Heaviside-style source profiles,
- optional smooth source turn-on/turn-off envelopes,
- optional ADE plasma response,
- optional PML absorbing boundaries,
- optional electric and magnetic divergence cleaning,
- optional normalized-unit operation,
- optional IMAS/IMASH loading of plasma profiles and magnetic fields,
- AMReX `MultiFab` storage and AMReX plotfile output,
- MPI domain decomposition through AMReX,
- optional CUDA compilation when AMReX and the project are built with CUDA.

---

## Code organization

Typical public headers are organized as follows:

```text
include/
├── fdtd/
│   ├── EHSolver.H
│   ├── EHJSolver.H
│   ├── PMLSolver.H
│   └── DivCleaning.H
├── io/
│   ├── Config.H
│   ├── ConfigHelpers.H
│   ├── FieldOutputManager.H
│   └── ImasReader.H
├── sources/
│   ├── Source.H
│   ├── PointSource.H
│   └── ADEPlasma.H
└── utils/
    ├── Constants.H
    ├── HelperFunctions.H
    └── Stencils.H
```

The main conceptual components are:

| Component | Purpose |
|---|---|
| `EHSolver` | Evolves the electromagnetic fields `E` and `H`. |
| `EHJSolver` | Extends the electromagnetic update with current-density coupling through `J`. |
| `PMLSolver` | Implements absorbing-layer damping near the physical boundaries. |
| `DivCleaning` | Applies divergence-cleaning corrections to reduce electric and/or magnetic divergence error. |
| `ADEPlasma` | Evolves an electron momentum response and accumulates plasma current into the electromagnetic solver. |
| `Source` | Abstract base class for time-dependent electromagnetic sources. |
| `PointSource` | Localized finite-width source deposited into `J`, `E`, or `H`. |
| `Config` | Reads and validates runtime options from the AMReX ParmParse input file. |
| `FieldOutputManager` | Registers fields and writes AMReX plotfiles. |
| `ImashSetup` / `ImasReader` | Build an FDTD setup from IMAS/IMASH plasma and magnetic-field data when enabled. |

---

## Theory and numerical model

### Maxwell system

The electromagnetic part of the code advances Maxwell's equations in the form

\f[
\mu_0\frac{\partial \mathbf{H}}{\partial t}
=
-\nabla \times \mathbf{E},
\f]

\f[
\epsilon_0\frac{\partial \mathbf{E}}{\partial t}
=
\nabla \times \mathbf{H}
,
\f]

or with the normalized equivalents when `init.normalize_units = 1`.

The solver stores and evolves the electromagnetic fields on AMReX `MultiFab` data structures.  The field update is FDTD-like: spatial derivatives are represented by finite-difference curl stencils, and the timestep is constrained by a CFL condition based on the cell spacing and the electromagnetic wave speed.

The timestep is controlled by

```text
init.cfl
```

and the simulation stops when either

```text
init.final_time
```

is reached or

```text
init.max_step
```

is reached.  A value of `init.max_step = -1` disables the step-count limit.

---

### Current coupling

When current-density coupling is active, Ampere's law contains the current source term

\f[
\epsilon_0\frac{\partial \mathbf{E}}{\partial t}
=
\nabla \times \mathbf{H}-\mathbf{J},
\f]
together with charge conservation
\f[
\partial_t\rho+\nabla\cdot\mathbf J=0.
\f]
The current field `J` can receive contributions from:

1. user-defined electromagnetic sources,
2. the ADE plasma model,
3. optional PML/plasma damping terms depending on the chosen configuration.

A source can target one of three fields:

```text
source.N.type = J
source.N.type = E
source.N.type = H
```

A `J` source is physically a current-density drive.  Direct `E` or `H` sources are useful for controlled numerical tests, wave launching, or diagnostics.

---

### ADE plasma response

The auxiliary differential equation (ADE) plasma model evolves an electron momentum-like response and accumulates a current into the electromagnetic solver.  Conceptually, the model follows a cold/warm electron-fluid response of the form

\f[
\frac{d\mathbf{p}_e}{dt}
=
q_e
\mathbf{E}
+
\mathbf{p}_e \times \mathbf{\Omega}_e
-
\nu \mathbf{p}_e,
\f]

where 

\f[
\mathbf {\Omega}_e=\frac{q_e\mathbf B}{m_e}
\f]

and a current contribution

\f[
\mathbf{J}_e
=
\frac{q_e n_e}{m_e} \mathbf{p}_e.
\f]
Note that the current density equation above is always true, even if all quantities vary in space and time. In the implementation, the ADE plasma object stores fields such as electron density `ne`, temperature `Te`, collision frequency, conductivity/activity mask, background magnetic field `B`, and electron momentum `p`.  The ADE model reads the electromagnetic field and deposits the resulting plasma current into the coupled `EHJSolver`.

For reflectometry in ITER, the electron mass is dependent on the plasma temperature as
\f[
    m_e=m_{0e}\sqrt{1+5\frac{eT_e}{m_{0e}c^2}},
\f]
where \f[T_e(t)\f] is the local electron temperature in eV. Thus the exact time advance can be written as
\f[
    \textbf p_e(t)
    =
    e^{\mathbb A t}\textbf p_0
    +
    q_e
    \left[
        c_0\mathbb I
        +
        c_1[\boldsymbol\Omega_e]_\times
        +
        c_2[\boldsymbol\Omega_e]_\times^2
    \right]
    \textbf E .
\f]

The coefficients are
\f[
    \left\{
    \begin{aligned}
    c_0
    &=
    \frac{1-e^{-\nu t}}{\nu},
    \\[0.5em]
    c_1
    &=
    \frac{
        \Omega_e
        -
        e^{-\nu t}
        \left(
            \nu\sin\theta+\Omega_e\cos\theta
        \right)
    }
    {
        \Omega_e
        \left(
            \nu^2+\Omega_e^2
        \right)
    },
    \\[0.5em]
    c_2
    &=
    \frac{1}{\Omega_e^2}
    \left[
        \frac{1-e^{-\nu t}}{\nu}
        -
        \frac{
            \nu
            +
            e^{-\nu t}
            \left(
                -\nu\cos\theta+\Omega_e\sin\theta
            \right)
        }
        {
            \nu^2+\Omega_e^2
        }
    \right].
    \end{aligned}
    \right.
\f]
with 
\f[
    e^{\mathbb A t}
    =
    e^{-\nu t}
    \left[
        \mathbb I
        +
        \frac{\sin\theta}{\Omega_e}
        [\boldsymbol\Omega_e]_\times
        +
        \frac{1-\cos\theta}{\Omega_e^2}
        [\boldsymbol\Omega_e]_\times^2
    \right],
\f]
where 
\f[
	\mathbb I=	
		\begin{bmatrix}
		1&0&0\\0&1&0\\0&0&1
		\end{bmatrix},
	\qquad
    [\boldsymbol\Omega_e]_\times
    =
    \begin{bmatrix}
        0&\Omega_{e_z}&-\Omega_{e_y}\\
        -\Omega_{e_z}&0&\Omega_{e_x}\\
        \Omega_{e_y}&-\Omega_{e_x}&0
    \end{bmatrix},
    \qquad
    \Omega_e=\sqrt{\boldsymbol\Omega_e\cdot\boldsymbol\Omega_e},
	\qquad
    \theta=\Omega_e t.
\f]

The plasma model can be initialized in two ways:

1. from manually specified/default plasma data,
2. from IMAS/IMASH-derived equilibrium and profile data.

---

### PML absorbing layer

The PML damps outgoing waves near selected domain boundaries.  It is enabled with

```text
init.use_pml = 1
```

and controlled by

```text
init.pml_thickness
init.pml_sigma
```

where `init.pml_thickness` is the layer thickness and `init.pml_sigma` controls the damping strength.

For plasma-coupled runs, the input may also include

```text
init.plasma_sigma
```

which scales the plasma conductivity mask associated with PML damping.

---

### Divergence cleaning

The code can apply divergence-cleaning corrections for electric and magnetic fields:

```text
init.use_e_cleaner = 1
init.use_h_cleaner = 0
```
By default, the magnetic field cleaner is turned off since the FDTD algorithm converses
\f[
\nabla\cdot\mathbf H=0
\f]
Electric field cleaning is useful when numerical current deposition, sources, or boundary treatments introduce error in the discrete Gauss-law constraint.  Magnetic cleaning can be enabled for tests where numerical magnetic divergence error needs active control.

---

### Normalized units

The code supports dimensional and normalized operation:

```text
init.normalize_units = 0
```

keeps the input and output in physical/code units, while

```text
init.normalize_units = 1
```

converts selected quantities using the normalization constants defined in `include/utils/Constants.H`.

Typical normalization constants include:

```text
B_0, L_0, n_0, E_0, H_0, t_0, J_0, rho_0
```

Use normalized mode consistently: source frequencies, source positions, IMAS-derived quantities, and output interpretation depend on whether normalization is enabled.

---

## Dependencies

Required:

- CMake 3.25 or newer,
- C++17 compiler,
- MPI,
- AMReX,
- Ninja or Make.

Optional:

- CUDA toolkit and compatible NVIDIA driver,
- Doxygen and Graphviz for documentation,
- IMASH/IMAS dependencies for ITER or IMAS equilibrium-based runs.

---

## Installing AMReX

Clone AMReX:

```bash
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex
```

Create a build directory:

```bash
mkdir build
cd build
```

### CPU/MPI AMReX build

```bash
cmake .. \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DAMReX_MPI=ON \
  -DAMReX_SPACEDIM=3 \
  -DCMAKE_INSTALL_PREFIX=$HOME/amrex_CPU/installdir

ninja
ninja install
```

### CUDA/MPI AMReX build

First identify the GPU compute capability:

```bash
nvidia-smi --query-gpu=compute_cap --format=csv
```

You can also inspect the architectures supported by `nvcc`:

```bash
nvcc --list-gpu-arch
```

Then configure AMReX with an explicit architecture.  For example, for compute capability 8.6:

```bash
cmake .. \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DAMReX_MPI=ON \
  -DAMReX_GPU_BACKEND=CUDA \
  -DAMReX_CUDA=ON \
  -DAMReX_SPACEDIM=3 \
  -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DCMAKE_INSTALL_PREFIX=$HOME/amrex_GPU/installdir

ninja
ninja install
```

Avoid relying on an automatic CUDA architecture setting if your project configuration requires the architecture to be stated explicitly.

---

## Building reflectomITER

Clone or enter the project directory:

```bash
cd reflectomITER
```

Create a clean build directory:

```bash
rm -rf build
mkdir build
```

### CPU build

Use this when CUDA is disabled:

```bash
cmake -S . -B build \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_CUDA=OFF \
  -DUSE_IMAS=OFF \
  -DAMReX_DIR=$HOME/amrex_CPU/installdir/lib/cmake/AMReX

cmake --build build
```
This line
```bash
  -DAMReX_DIR=$HOME/amrex_CPU/installdir/lib/cmake/AMReX
```
can go if AMReX installed globaly or as a module.
### CUDA build

Use this when AMReX was built with CUDA:

```bash
cmake -S . -B build \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_CUDA=ON \
  -DUSE_IMAS=OFF \
  -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DAMReX_DIR=$HOME/amrex_GPU/installdir/lib/cmake/AMReX

cmake --build build
```
This line
```bash
  -DAMReX_DIR=$HOME/amrex_GPU/installdir/lib/cmake/AMReX
```
can go if AMReX installed globaly or as a module.

### IMAS/IMASH build

If the project is configured to fetch IMASH through CMake, enable IMAS support with:

```bash
cmake -S . -B build \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_CUDA=ON \
  -DUSE_IMAS=ON \
  -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DAMReX_DIR=$HOME/amrex_GPU/installdir/lib/cmake/AMReX

cmake --build build
```

If you are on an HPC system, load the required compiler, MPI, CUDA, AMReX, IMAS, and IMASH modules before configuring the project.

### iWrap wrapper build and usage

The `iwrap/` directory contains an actor-style wrapper scaffold for launching `reflectomITER` from an iWrap-oriented workflow.

There are currently two wrapper paths:

```text
Python wrapper
  recommended for local WSL/debugging
  does not require the C++ MUSCLE3/libmuscle development stack

C++ iWrap adapter
  useful for testing iWrap-style compiled entry points
  closer to the future compiled/library actor path
```

The production ITER target is expected to move toward a library actor, but the executable wrapper remains useful for validating input generation, source setup, IMAS configuration, and MPI launch behavior.

#### Python wrapper prerequisites

The Python wrapper only requires a working Python interpreter and a compiled `reflectomITER` executable.

Optional but recommended:

```bash
python -m pip install --upgrade pip
```

The Python wrapper can be run without installing iWrap itself.

#### Python wrapper smoke test

From the repository root, build `reflectomITER` first using one of the build modes above. Then run:

```bash
cd iwrap

python3 run_reflectomITER_actor.py
```

The run configuration is currently programmed in:

```text
iwrap/run_reflectomITER_actor.py
```

The reusable wrapper package lives in:

```text
iwrap/reflectomITER_actor/
```

The Python wrapper generates an input file from a fully populated template by replacing selected parameters such as:

```text
init.output_dir
init.max_step
init.final_time
init.cfl
init.number_of_outputs
init.nsources
source.N.frequency
source.N.amplitude
```

The default generated input and report are written to paths configured in `run_reflectomITER_actor.py`.

#### Python wrapper structure

```text
iwrap/
├── run_reflectomITER_actor.py
├── reflectomITER_python_actor.py
├── reflectomITER_iwrap.yml
├── test_iwrap_adapter.cpp
└── reflectomITER_actor/
    ├── __init__.py
    ├── actor.py
    ├── config.py
    ├── input_generation.py
    ├── launcher.py
    └── sources.py
```

The file `run_reflectomITER_actor.py` is the high-level run script. It defines the base run parameters, source descriptions, generated input path, launcher, and executable path.

The package `reflectomITER_actor/` contains reusable code for:

```text
input-file generation
source-to-parameter conversion
MPI/subprocess launch
actor-style lifecycle functions
```

#### Line endings on WSL

If a Python script fails with:

```text
/usr/bin/env: ‘python3\r’: No such file or directory
```

then the file has Windows CRLF line endings. Fix all Python files under `iwrap/` with:

```bash
cd iwrap
find . -name "*.py" -print0 | xargs -0 sed -i 's/\r$//'
chmod +x run_reflectomITER_actor.py reflectomITER_python_actor.py
```

To prevent this from recurring, keep a repository-level `.gitattributes` file with:

```text
*.py text eol=lf
*.sh text eol=lf
*.yml text eol=lf
*.yaml text eol=lf
*.cpp text eol=lf
*.H text eol=lf
*.h text eol=lf
*.md text eol=lf
```

#### C++ iWrap adapter smoke test

The C++ adapter exposes iWrap-style lifecycle functions from:

```text
include/io/IWrapActor.H
src/io/IWrapActor.cpp
```

A simple local smoke-test executable can be built manually with:

```bash
g++ -std=c++17 -Iinclude \
    src/io/IWrapActor.cpp \
    iwrap/test_iwrap_adapter.cpp \
    -o iwrap/test_iwrap_adapter
```

Then run:

```bash
cd iwrap

REFLECTOMITER_LAUNCHER=mpirun \
REFLECTOMITER_LAUNCHER_ARGS="-np 2" \
REFLECTOMITER_EXECUTABLE=../build/reflectomITER \
REFLECTOMITER_INPUT=../inputs/input_IMAS.txt \
./test_iwrap_adapter
```

Adjust `REFLECTOMITER_EXECUTABLE` if the solver executable was built in another directory, for example:

```bash
REFLECTOMITER_EXECUTABLE=../build_iwrap/reflectomITER
```

#### Installing iWrap locally

For local development, clone and install iWrap from source:

```bash
cd ~
git clone https://github.com/iterorganization/iWrap.git
cd iWrap

python -m pip install --upgrade pip setuptools wheel setuptools_scm
python -m pip install -v -e .
```

Check that the command is available:

```bash
iwrap --list-actor-types
```

If `pip install iwrap` does not work, install from the GitHub source tree as shown above.

#### IMAS Access Layer visibility

The compiled iWrap path may require IMAS Access Layer packages visible through `pkg-config`.

For a local IMAS-Cpp installation, expose `al-cpp.pc` with:

```bash
export PKG_CONFIG_PATH="$HOME/IMAS-Cpp/installdir/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
```

Verify:

```bash
pkg-config --modversion al-cpp
pkg-config --cflags al-cpp
pkg-config --libs al-cpp
```

On an ITER or HPC system, this should normally be handled by the official IMAS/iWrap environment setup or modules.

#### C++ MUSCLE3/libmuscle requirements

Generating a full compiled C++ iWrap/MUSCLE3 actor may additionally require:

```bash
pkg-config --modversion ymmsl
pkg-config --modversion libmuscle_mpi
```

and the C++ header:

```cpp
#include <libmuscle/libmuscle.hpp>
```

If these are not available on WSL, use the Python wrapper path locally and build the full C++ actor on an ITER/Gateway environment where IMAS, iWrap, MUSCLE3, MPI, and compiler modules are provided.

#### Generating a Python iWrap actor

Once iWrap is installed, the Python actor scaffold can be generated with:

```bash
cd iwrap

iwrap --actor-type python \
  -a reflectomITER_actor \
  -f reflectomITER_iwrap.yml
```

The exact actor type names available in a given iWrap installation can be checked with:

```bash
iwrap --list-actor-types
```

If local iWrap template generation fails because the Access Layer version is undefined, make sure the IMAS environment is sourced. On WSL, this may require patching local iWrap templates or using the Python wrapper directly without generating the actor.

#### Current recommended workflow

For local WSL development:

```text
1. Build reflectomITER normally.
2. Run iwrap/run_reflectomITER_actor.py.
3. Validate generated input files and reports.
4. Validate MPI launch behavior.
5. Commit wrapper/input-generation changes.
```

For ITER production integration:

```text
1. Use the executable wrapper to validate physics setup and input generation.
2. Refactor reflectomITER into a callable core library.
3. Expose the core through an iWrap library actor.
4. Add IDS input/output handling at the actor boundary.
```


---

## Running the code

The executable name depends on the `add_executable(...)` line in `CMakeLists.txt`.  In the current CMake setup, the target is commonly

```text
reflectomITER
```

so the executable is typically:

```bash
./build/reflectomITER
```

If your local branch still names the executable `fdtd`, replace `reflectomITER` by `fdtd` in the examples below.

### Serial run

```bash
./build/reflectomITER input.txt
```

### MPI run

```bash
mpirun -np 4 ./build/reflectomITER input.txt
```

or, from inside the build directory:

```bash
cd build
mpirun -np 4 ./reflectomITER ../input.txt
```

For an IMAS-backed run:

```bash
mpirun -np 4 ./build/reflectomITER input_IMAS.txt
```

Note that the number of MPI tasks must match the number of GPU you have on the machine if you are running with GPUs.

---

## Input-file structure

The code uses AMReX `ParmParse`-style input files.  Parameters are grouped by prefixes:

| Prefix | Meaning |
|---|---|
| `init.*` | Runtime, solver, PML, output, plasma, and cleaning options. |
| `input.*` | Manual-domain geometry and grid options. |
| `imas.*` | IMAS/IMASH equilibrium-loading and generated-domain options. |
| `source.N.*` | Parameters for source number `N`. |

Important common options:

| Option | Meaning |
|---|---|
| `init.nghost` | Number of ghost cells. |
| `init.normalize_units` | Use normalized units if nonzero. |
| `init.cfl` | CFL number used for timestep selection. |
| `init.final_time` | Final simulation time. |
| `init.max_step` | Maximum number of steps; `-1` disables the limit. |
| `init.number_of_outputs` | Number of requested output intervals. |
| `init.use_plasma` | Enable ADE plasma response. |
| `init.use_pml` | Enable PML damping. |
| `init.use_e_cleaner` | Enable electric divergence cleaning. |
| `init.use_h_cleaner` | Enable magnetic divergence cleaning. |
| `init.output_dir` | Root directory for output files. |
| `init.nsources` | Number of source blocks to read. |

Source options:

| Option | Meaning |
|---|---|
| `source.N.enabled` | Enable or disable source `N`. |
| `source.N.type` | Target field: `J`, `E`, or `H`. |
| `source.N.position` | Source center in physical/code coordinates. |
| `source.N.direction` | Source polarization vector. |
| `source.N.shape` | Anisotropic profile stretch factors. |
| `source.N.amplitude` | Source amplitude. |
| `source.N.frequency` | Oscillation frequency. |
| `source.N.phase` | Phase in degrees. |
| `source.N.width` | Spatial width of the source. |
| `source.N.profile` | `gaussian`, `supergaussian`, or `heaviside`. |
| `source.N.supergaussian_order` | Order for super-Gaussian profile. |
| `source.N.t_on` | Optional smooth turn-on time. |
| `source.N.t_off` | Optional smooth turn-off time. |
| `source.N.t_rise` | Optional tanh rise/fall time. |

---

## Example 1: manual-domain run

The manual-domain input file disables IMAS and specifies the computational domain directly.

Example command:

```bash
mpirun -np 1 ./build/reflectomITER input.txt
```

Key features of this example:

- dimensional/code-unit run,
- 2-D-like domain with one periodic `y` cell,
- PML enabled,
- ADE plasma enabled,
- electric divergence cleaning enabled,
- two requested sources, with the first source depositing into `E`.

Representative input:

```text
# ============================================================
# Common runtime / solver
# ============================================================
init.nghost = 1
init.normalize_units = 0

init.cfl = 0.9
init.final_time = 2.5
init.max_step = -1
init.number_of_outputs = 10

init.use_plasma = 1
init.use_e_cleaner = 1
init.use_h_cleaner = 0

init.output_dir = "/home/user/data"

# ============================================================
# PML Setup
# ============================================================
init.use_pml = 1
init.pml_thickness = 30
init.pml_sigma = 50

# ============================================================
# Manual input domain
# ============================================================
input.scale = 10
input.n_cell = 32 1 32

input.prob_lo = -1.0 -1e-2 -1.0
input.prob_hi =  1.0  1e-2  1.0

input.is_periodic = 0 1 0

# ============================================================
# IMAS disabled
# ============================================================
imas.use_imas = 0

# ============================================================
# Sources
# ============================================================
init.nsources = 2

source.1.enabled = 1
source.1.type = E
source.1.position = 0.75 0.0 0.0
source.1.direction = 0.0 1.0 0.0
source.1.shape = 1.0 1.0 1.0
source.1.amplitude = 1e-1
source.1.frequency = 10
source.1.phase = 0.0
source.1.width = 0.02
source.1.profile = supergaussian
source.1.supergaussian_order = 2.0
```

The cell count is read from

```text
input.n_cell = 32 1 32
```

and then scaled by

```text
input.scale = 10
```

except for collapsed directions with one cell.  Therefore this input represents a high-resolution `x-z` run with a single periodic cell in `y`.
For three dimensional runs replace the domain block with its three dimensional version
```text
# ============================================================
# Manual input domain
# ============================================================
input.scale = 10
input.n_cell = 32 32 32

input.prob_lo = -1.0 -1.0 -1.0
input.prob_hi =  1.0  1.0  1.0

input.is_periodic = 0 0 0
```

---

## Example 2: IMAS/IMASH run

The IMAS/IMASH input file builds the simulation geometry and plasma/background-field setup from an IMAS data source.

Example command:

```bash
mpirun -np 4 ./build/reflectomITER input_IMAS.txt
```

Key features of this example:

- IMAS loading enabled,
- normalized units enabled,
- core profiles, edge profiles, and magnetic field loaded,
- domain generated from IMAS data and a toroidal wedge angle,
- PML enabled with stronger damping,
- ADE plasma enabled,
- current-density source at GHz frequency,
- smooth source turn-on/off envelope.

Representative input:

```text
# ============================================================
# IMASH / equilibrium loading
# ============================================================
imas.use_imas = 1
imas.imas_path = /home/iter/134174/117
imas.imas_time = 60.0

imas.load_core = 1
imas.load_bfield = 1
imas.load_edge_on_grid = 1
imas.load_wall = 0
imas.smooth_rho = 0.01
imas.smooth_theta = 0.2

imas.merge_density = 1
imas.merge_temperature = 1

# ============================================================
# Domain generation from IMASH data
# ============================================================
imas.toroidal_phi_deg = 3
imas.points_per_wavelength = 18.0
imas.pad_R = 0.0
imas.pad_Z = 0.0
imas.reference_frequency = -1.0

imas.R_min = 6.0
imas.R_max = 8.6
imas.Z_min = 0.0
imas.Z_max = 1.0

# ============================================================
# Grid / solver
# ============================================================
init.nghost = 1
init.normalize_units = 1

# ============================================================
# PML Setup
# ============================================================
init.pml_thickness = 80
init.pml_sigma = 1500
init.plasma_sigma = 200
init.use_pml = 1

# ============================================================
# Time stepping
# ============================================================
init.cfl = 0.45
init.final_time = 15e-9
init.max_step = -1
init.number_of_outputs = 15

# ============================================================
# Feature switches
# ============================================================
init.use_plasma = 1
init.use_e_cleaner = 1
init.use_h_cleaner = 0

# ============================================================
# Output
# ============================================================
init.output_dir = "/home/user/data"

# ============================================================
# Sources
# ============================================================
init.nsources = 1

source.1.enabled = 1
source.1.type = J
source.1.position = 8.5 0.0 0.5
source.1.direction = 0.0 1.0 0.0
source.1.shape = 0.05 1.0 1.0
source.1.amplitude = 1e-3
source.1.frequency = 37e9
source.1.phase = 0.0
source.1.width = .05
source.1.profile = gaussian
source.1.supergaussian_order = 4.0
source.1.t_on = 0.0
source.1.t_off = 33e-9
source.1.t_rise = .1e-9
```

The setting

```text
imas.reference_frequency = -1.0
```

means that the reference frequency used for grid generation should be inferred from the enabled source frequencies.  In this example, the enabled source frequency is

```text
source.1.frequency = 37e9
```

and the requested resolution is

```text
imas.points_per_wavelength = 18.0
```

The settings

```text
imas.R_min = 6.0
imas.R_max = 8.6
imas.Z_min = 0.0
imas.Z_max = 1.0
```

restrict the generated IMAS-derived computational domain.

---

## Output

The code writes AMReX plotfiles under

```text
init.output_dir
```

For example:

```text
init.output_dir = "/home/user/data"
```

`FieldOutputManager` registers solver and diagnostic fields, applies optional unit scaling, and writes AMReX-compatible output.  The output can typically be inspected using AMReX-aware visualization workflows, VisIt, ParaView with appropriate readers, or custom Python post-processing.

The number of output intervals is controlled by

```text
init.number_of_outputs
```

The implementation may also write metadata files that help visualization tools load a time series of plotfiles.

---

## Generating Doxygen documentation

Install Doxygen and Graphviz:

```bash
sudo apt install doxygen graphviz
```

or on an HPC system:

```bash
module avail Doxygen
module load Doxygen
```

Configure with documentation enabled:

```bash
cmake -S . -B build -DBUILD_DOCS=ON
```

Generate the documentation:

```bash
cmake --build build --target docs
```

The HTML entry point is usually:

```text
build/docs/html/index.html
```

Open it locally with:

```bash
firefox build/docs/html/index.html
```

A useful CMake Doxygen block is:

```cmake
option(BUILD_DOCS "Build Doxygen documentation" ON)

if(BUILD_DOCS)
    find_package(Doxygen)

    if(DOXYGEN_FOUND)
        set(DOXYGEN_PROJECT_NAME "reflectomITER")
        set(DOXYGEN_PROJECT_BRIEF "AMReX-based FDTD Maxwell solver")
        set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/docs")

        set(DOXYGEN_GENERATE_HTML YES)
        set(DOXYGEN_GENERATE_LATEX NO)
        set(DOXYGEN_GENERATE_TREEVIEW YES)

        set(DOXYGEN_EXTRACT_ALL YES)
        set(DOXYGEN_EXTRACT_CLASSES YES)
        set(DOXYGEN_EXTRACT_PRIVATE YES)
        set(DOXYGEN_EXTRACT_STATIC YES)
        set(DOXYGEN_EXTRACT_LOCAL_CLASSES YES)
        set(DOXYGEN_HIDE_UNDOC_CLASSES NO)
        set(DOXYGEN_HIDE_UNDOC_MEMBERS NO)

        set(DOXYGEN_RECURSIVE YES)

        set(DOXYGEN_FILE_PATTERNS
            *.md
            *.h
            *.H
            *.hh
            *.hpp
            *.hxx
            *.c
            *.cc
            *.cpp
            *.cxx
            *.cu
        )

        set(DOXYGEN_ENABLE_PREPROCESSING YES)
        set(DOXYGEN_MACRO_EXPANSION YES)
        set(DOXYGEN_EXPAND_ONLY_PREDEF NO)

        set(DOXYGEN_PREDEFINED
            AMREX_GPU_DEVICE=
            AMREX_GPU_HOST_DEVICE=
            AMREX_FORCE_INLINE=inline
            AMREX_SPACEDIM=3
        )

        if(EXISTS "${CMAKE_SOURCE_DIR}/README.md")
            set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_SOURCE_DIR}/README.md")
        endif()

        doxygen_add_docs(
            docs
            "${CMAKE_SOURCE_DIR}/README.md"
            "${CMAKE_SOURCE_DIR}/include"
            "${CMAKE_SOURCE_DIR}/src"
            COMMENT "Generating Doxygen documentation"
        )
    else()
        message(WARNING "Doxygen not found: documentation target 'docs' will not be created")
    endif()
endif()
```

---

## Troubleshooting

### `ninja: error: unknown target 'docs'`

The Doxygen target was not created.  Reconfigure after installing Doxygen:

```bash
cmake -S . -B build -DBUILD_DOCS=ON
cmake --build build --target help | grep -i docs
```

If the target still does not appear, check that `find_package(Doxygen)` succeeds during configuration.

---

### Doxygen runs but no classes are shown

Make sure these options are set:

```cmake
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_EXTRACT_CLASSES YES)
set(DOXYGEN_HIDE_UNDOC_CLASSES NO)
set(DOXYGEN_FILE_PATTERNS *.H *.h *.hpp *.cpp *.cu)
set(DOXYGEN_RECURSIVE YES)
```

For AMReX/CUDA macros, also define:

```cmake
set(DOXYGEN_PREDEFINED
    AMREX_GPU_DEVICE=
    AMREX_GPU_HOST_DEVICE=
    AMREX_FORCE_INLINE=inline
    AMREX_SPACEDIM=3
)
```

---

### CUDA is still detected when `USE_CUDA=OFF`

Do not list CUDA unconditionally in the `project()` language list.  Put the `USE_CUDA` option before `project()`:

```cmake
cmake_minimum_required(VERSION 3.25)

option(USE_CUDA "Enable CUDA GPU support" OFF)

if(USE_CUDA)
    project(reflectomITER LANGUAGES C CXX CUDA)
else()
    project(reflectomITER LANGUAGES C CXX)
endif()
```

Then configure from a clean build directory:

```bash
rm -rf build
cmake -S . -B build -DUSE_CUDA=OFF
```

---

### CUDA driver/runtime mismatch

If you see an error such as

```text
CUDA driver version is insufficient for CUDA runtime version
```

then the runtime used to build AMReX or the application is newer than the NVIDIA driver available on the node.  Check the node driver with:

```bash
nvidia-smi
```

On a cluster, check the driver on the compute node, not only the login node.

---

### AMReX compiler mismatch

If AMReX was built with one compiler family, build this code with a compatible compiler family.  For example, if AMReX was built with Intel compilers, configure the application with the same compilers:

```bash
export CC=icc
export CXX=icpc
cmake -S . -B build ...
```

or use the equivalent modern Intel compiler names if your environment uses `icx`/`icpx`.


