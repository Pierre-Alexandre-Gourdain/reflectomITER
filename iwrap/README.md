# reflectomITER iWrap wrapper scaffold

This directory contains the initial wrapper scaffold for running `reflectomITER` from an iWrap/actor-style workflow.

The current goal is simple:

```text
wrapper
  -> optional launcher: mpirun, srun, ...
  -> reflectomITER executable
  -> existing AMReX ParmParse input file
```

The wrapper does not yet generate the input file from IDS data. It currently launches `reflectomITER` using an existing input file.

## Files

```text
iwrap/reflectomITER_python_actor.py   Python launch wrapper for WSL/debugging
iwrap/test_iwrap_adapter.cpp          C++ adapter smoke-test driver
iwrap/reflectomITER_iwrap.yml         Initial iWrap code-description scaffold
include/io/IWrapActor.H               C++ iWrap-style adapter declarations
src/io/IWrapActor.cpp                 C++ iWrap-style adapter implementation
```

## Python wrapper

The Python wrapper reads these environment variables:

```text
REFLECTOMITER_EXECUTABLE       path to reflectomITER
REFLECTOMITER_INPUT            path to input file
REFLECTOMITER_LAUNCHER         optional launcher, e.g. mpirun
REFLECTOMITER_LAUNCHER_ARGS    launcher arguments, e.g. "-np 2"
REFLECTOMITER_WORKDIR          working directory
REFLECTOMITER_IWRAP_LOG        log file
REFLECTOMITER_DRY_RUN          if 1, print command without running
```

Example dry run:

```bash
cd ~/reflectomITER/iwrap

REFLECTOMITER_DRY_RUN=1 \
REFLECTOMITER_LAUNCHER=mpirun \
REFLECTOMITER_LAUNCHER_ARGS="-np 2" \
REFLECTOMITER_EXECUTABLE=../build_iwrap/reflectomITER \
REFLECTOMITER_INPUT=../inputs/input.txt \
./reflectomITER_python_actor.py
```

Example real run:

```bash
REFLECTOMITER_LAUNCHER=mpirun \
REFLECTOMITER_LAUNCHER_ARGS="-np 2" \
REFLECTOMITER_EXECUTABLE=../build_iwrap/reflectomITER \
REFLECTOMITER_INPUT=../inputs/input.txt \
./reflectomITER_python_actor.py
```

Use `../build/reflectomITER` instead if that is where the solver executable was built.

The wrapper writes its log to:

```text
reflectomITER_iwrap_python.log
```

## C++ adapter smoke test

Compile from the repository root:

```bash
g++ -std=c++17 -Iinclude \
    src/io/IWrapActor.cpp \
    iwrap/test_iwrap_adapter.cpp \
    -o iwrap/test_iwrap_adapter
```

Run:

```bash
cd iwrap

REFLECTOMITER_LAUNCHER=mpirun \
REFLECTOMITER_LAUNCHER_ARGS="-np 2" \
REFLECTOMITER_EXECUTABLE=../build_iwrap/reflectomITER \
REFLECTOMITER_INPUT=../inputs/input.txt \
./test_iwrap_adapter
```

## iWrap notes

The Python wrapper is the easiest local WSL path.

The compiled C++ iWrap path may require IMAS/MUSCLE3 packages visible through `pkg-config`, for example:

```bash
pkg-config --modversion al-cpp
pkg-config --modversion ymmsl
pkg-config --modversion libmuscle_mpi
```

For local testing, `al-cpp` was exposed with:

```bash
export PKG_CONFIG_PATH="/home/pag/IMAS-Cpp/installdir/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
```

Full actor generation is expected to be easier on the ITER/Gateway environment after sourcing the official IMAS/iWrap setup scripts.

## Current status

Validated locally:

```text
reflectomITER_python_actor.py
  -> mpirun -np 2
  -> reflectomITER
  -> existing input file
```

### iWrap wrapper build and usage

The `iwrap/` directory contains an actor-style wrapper scaffold for launching `reflectomITER` from an iWrap-oriented workflow.

There are currently two wrapper paths:

```text
Python wrapper
  recommended for local WSL/debugging
  does not require the C++ MUSCLE3/libmuscle development stack

C++ iWrap adapter
  exposes iWrap-style C++ lifecycle functions
  useful for testing compiled actor entry points
  closer to the future compiled/library actor path
```

The production ITER target is expected to move toward a library actor, but the executable wrapper remains useful for validating input generation, source setup, IMAS configuration, and MPI launch behavior.

---

#### Python wrapper prerequisites

The Python wrapper only requires a working Python interpreter and a compiled `reflectomITER` executable.

Optional but recommended:

```bash
python -m pip install --upgrade pip
```

The Python wrapper can be run without installing iWrap itself.

---

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

---

#### C++ iWrap adapter files

The C++ adapter lives in:

```text
include/io/IWrapActor.H
src/io/IWrapActor.cpp
```

It exposes iWrap-style lifecycle functions:

```cpp
init_code(...)
code_step(...)
clean_up(...)
get_code_state(...)
restore_code_state(...)
get_timestamp_cpp(...)
```

The implementation is currently an executable launcher. It builds a command such as:

```text
mpirun -np 2 ../build/reflectomITER ../inputs/input_IMAS.txt
```

and runs it from the wrapper process.

The C++ adapter is independent of the full solver link. It does not directly link AMReX or the `reflectomITER` solver core at this stage; it only launches the already-built executable.

---

#### Compile the C++ iWrap adapter smoke test independently

From the repository root:

```bash
g++ -std=c++17 -Iinclude \
    src/io/IWrapActor.cpp \
    iwrap/test_iwrap_adapter.cpp \
    -o iwrap/test_iwrap_adapter
```

This builds only the small adapter test executable:

```text
iwrap/test_iwrap_adapter
```

It does not build the full solver. The full solver must already exist, for example:

```text
build/reflectomITER
```

or:

```text
build_iwrap/reflectomITER
```

Run a dry test:

```bash
cd iwrap

REFLECTOMITER_DRY_RUN=1 \
REFLECTOMITER_LAUNCHER=mpirun \
REFLECTOMITER_LAUNCHER_ARGS="-np 2" \
REFLECTOMITER_EXECUTABLE=../build/reflectomITER \
REFLECTOMITER_INPUT=../inputs/input_IMAS.txt \
./test_iwrap_adapter
```

Run for real:

```bash
REFLECTOMITER_LAUNCHER=mpirun \
REFLECTOMITER_LAUNCHER_ARGS="-np 2" \
REFLECTOMITER_EXECUTABLE=../build/reflectomITER \
REFLECTOMITER_INPUT=../inputs/input_IMAS.txt \
./test_iwrap_adapter
```

If the solver executable was built elsewhere, change the path:

```bash
REFLECTOMITER_EXECUTABLE=../build_iwrap/reflectomITER
```

The adapter log is written by default to:

```text
reflectomITER_iwrap.log
```

---

#### Compile the C++ iWrap adapter as a static library

iWrap's compiled wrapper path expects a native code library. The C++ adapter can be compiled into a static archive independently with:

```bash
cd iwrap

g++ -std=c++17 -fPIC -I../include \
    -c ../src/io/IWrapActor.cpp \
    -o IWrapActor.o

ar -cr libreflectomITER_iwrap_adapter.a IWrapActor.o
```

This produces:

```text
iwrap/libreflectomITER_iwrap_adapter.a
```

This archive contains only the adapter functions. It does not contain the full `reflectomITER` solver.

The corresponding iWrap YAML should point to the archive and header, for example:

```yaml
code_description:
  implementation:
    subroutines:
      init: init_code
      main:
        name: code_step
        arguments: []
      finalize: clean_up
      get_timestamp: get_timestamp_cpp
    programming_language: cpp
    data_dictionary_compliant: 3.39.0
    data_type: legacy
    code_path: ./libreflectomITER_iwrap_adapter.a
    include_path: ../include/io/IWrapActor.H

  documentation: "Executable-launch iWrap adapter for reflectomITER"

  settings:
    compiler_cmd: g++
    mpi_compiler_cmd: mpic++
    compiler_flags: -pthread -std=c++17 -I../include
```

Because iWrap's generated wrapper calls lifecycle functions in the global namespace, `include/io/IWrapActor.H` must expose global aliases for the namespaced implementation, for example:

```cpp
using reflectomiter::io::iwrap::init_code;
using reflectomiter::io::iwrap::clean_up;
using reflectomiter::io::iwrap::code_step;
using reflectomiter::io::iwrap::get_code_state;
using reflectomiter::io::iwrap::restore_code_state;
using reflectomiter::io::iwrap::get_timestamp_cpp;
```

---

#### Optional CMake target for the C++ adapter test

The adapter smoke test can also be added as an optional CMake target:

```cmake
option(BUILD_IWRAP_ADAPTER "Build reflectomITER iWrap adapter test driver" OFF)

if(BUILD_IWRAP_ADAPTER)
    add_library(reflectomITER_iwrap_adapter STATIC
        src/io/IWrapActor.cpp
    )

    target_include_directories(reflectomITER_iwrap_adapter PUBLIC include)

    add_executable(reflectomITER_iwrap_test
        iwrap/test_iwrap_adapter.cpp
    )

    target_link_libraries(reflectomITER_iwrap_test PRIVATE
        reflectomITER_iwrap_adapter
    )
endif()
```

Configure with:

```bash
cmake -S . -B build_iwrap_adapter \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_IWRAP_ADAPTER=ON

cmake --build build_iwrap_adapter
```

This builds the adapter test target, not necessarily the full solver unless the main project targets are also enabled in the same configuration.

---

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

---

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

---

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

---

#### Generate a Python iWrap actor

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

---

#### Generate a compiled C++ iWrap actor from the adapter archive

After creating:

```text
iwrap/libreflectomITER_iwrap_adapter.a
```

generate the actor with:

```bash
cd iwrap

iwrap --actor-type python \
  -a reflectomITER_actor \
  -f reflectomITER_iwrap.yml
```

or, for the MUSCLE3 C++ actor type if available and the required C++ dependencies are installed:

```bash
iwrap --actor-type muscle3-cpp \
  -a reflectomITER_actor \
  -f reflectomITER_iwrap.yml
```

The generated actor is installed by default under:

```text
$HOME/IWRAP_ACTORS/
```

Note that this compiled actor still launches the external `reflectomITER` executable. It is not yet a true linked solver-library actor.

---

#### Current recommended workflow

For local WSL development:

```text
1. Build reflectomITER normally.
2. Run iwrap/run_reflectomITER_actor.py.
3. Validate generated input files and reports.
4. Validate MPI launch behavior.
5. Optionally compile iwrap/test_iwrap_adapter.
6. Optionally build libreflectomITER_iwrap_adapter.a.
7. Commit wrapper/input-generation changes.
```

For ITER production integration:

```text
1. Use the executable wrapper to validate physics setup and input generation.
2. Refactor reflectomITER into a callable core library.
3. Expose the core through an iWrap library actor.
4. Add IDS input/output handling at the actor boundary.
```
