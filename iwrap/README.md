# reflectomITER iWrap wrapper build and usage

The `iwrap/` directory contains the current actor-wrapper scaffold for running `reflectomITER` from an iWrap-oriented workflow.

The current recommended structure is:

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

The current workflow is:

```text
run_reflectomITER_actor.py
  -> builds an ActorConfig
  -> defines base run parameters
  -> defines one or more source descriptions
  -> generates a reflectomITER input file
  -> launches reflectomITER through mpirun/srun or directly
```

The reusable code lives in the `reflectomITER_actor/` Python package. The run-specific physics setup should stay in `run_reflectomITER_actor.py`.

---

## Python actor wrapper

The Python wrapper is the recommended local development path. It does not require the full C++ MUSCLE3/libmuscle development stack.

From the repository root, build `reflectomITER` first using one of the normal build modes. Then run:

```bash
cd iwrap
python3 run_reflectomITER_actor.py
```

The run is configured directly in:

```text
iwrap/run_reflectomITER_actor.py
```

A typical run script creates an `ActorConfig`:

```python
config = ActorConfig(
    launcher="mpirun",
    launcher_arguments="-np 1",
    executable="../build/reflectomITER",

    input_template=Path("../inputs/input_IMAS.txt"),
    generated_input=Path("../inputs/input_IMAS_gen.txt"),

    working_directory=Path("."),
    log_file=Path("reflectomITER_iwrap_python.log"),
    dry_run=False,

    overrides={
        "init.output_dir": "../run_iwrap_test",
        "init.max_step": "100",
        "init.final_time": "15e-9",
        "init.cfl": "0.45",
        "init.number_of_outputs": "5",
    },
)
```

For a dry run, set:

```python
dry_run=True
```

The Python actor generates a new input file from a fully populated template by replacing selected parameters. The template can be a complete working `reflectomITER` input file.

---

## Multi-source setup

Multiple sources should be described in `run_reflectomITER_actor.py`.

The base source description can be written once:

```python
BASE_SOURCE = {
    "enabled": "1",
    "type": "J",
    "position": "8.5 0.0 0.5",
    "direction": "0.0 1.0 0.0",
    "shape": "0.05 1.0 1.0",
    "amplitude": "1e-3",
    "frequency": "3.7e9",
    "phase": "0.0",
    "width": ".05",
    "profile": "gaussian",
    "supergaussian_order": "4.0",
    "t_on": "0.0",
    "t_off": "33e-9",
    "t_rise": ".1e-9",
}
```

Then individual sources can be created by copying the base source and applying changes:

```python
SOURCES = [
    make_source(
        index=1,
        base_source=BASE_SOURCE,
        updates={
            "position": "8.5 0.0 0.5",
            "frequency": "35e9",
            "amplitude": "1e-3",
        },
    ),
    make_source(
        index=2,
        base_source=BASE_SOURCE,
        updates={
            "position": "8.2 0.0 0.5",
            "frequency": "37e9",
            "amplitude": "5e-4",
        },
    ),
    make_source(
        index=3,
        base_source=BASE_SOURCE,
        updates={
            "position": "8.0 0.0 0.5",
            "frequency": "40e9",
            "amplitude": "5e-4",
        },
    ),
]
```

The helper:

```python
sources_to_overrides(SOURCES)
```

generates parameters such as:

```text
init.nsources = 3
source.1.frequency = ...
source.2.frequency = ...
source.3.frequency = ...
```

These source parameters are merged with the base run overrides:

```python
overrides = {
    **BASE_RUN_OVERRIDES,
    **sources_to_overrides(SOURCES),
}
```

---

## Generated input files

The generated input file is written to the path specified by:

```python
generated_input=Path("../inputs/input_IMAS_gen.txt")
```

The input-generation module also writes a sidecar report:

```text
../inputs/input_IMAS_gen.txt.report
```

The report lists which parameters were replaced and which were appended.

Inspect the generated file with:

```bash
grep -E "init.nsources|source.[123].|init.output_dir|init.max_step|init.final_time|init.cfl" \
  ../inputs/input_IMAS_gen.txt

cat ../inputs/input_IMAS_gen.txt.report
```

---

## C++ iWrap adapter

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

At this stage, the C++ adapter launches an already-built `reflectomITER` executable. It does not link the full solver into the actor.

Compile the C++ smoke test from the repository root:

```bash
g++ -std=c++17 -Iinclude \
    src/io/IWrapActor.cpp \
    iwrap/test_iwrap_adapter.cpp \
    -o iwrap/test_iwrap_adapter
```

Run:

```bash
cd iwrap
./test_iwrap_adapter
```

This validates the C++ lifecycle functions independently of the full iWrap actor generator.

---

## Build the C++ adapter as a static archive

To use the compiled adapter with iWrap, build it as a static library:

```bash
cd iwrap

g++ -std=c++17 -fPIC -I../include \
    -c ../src/io/IWrapActor.cpp \
    -o IWrapActor.o

ar -cr libreflectomITER_iwrap_adapter.a IWrapActor.o
```

This creates:

```text
iwrap/libreflectomITER_iwrap_adapter.a
```

The archive contains only the adapter entry points. It does not contain the full `reflectomITER` solver.

Because iWrap's generated wrapper calls lifecycle functions in the global namespace, `include/io/IWrapActor.H` should expose global aliases for the namespaced implementation:

```cpp
using reflectomiter::io::iwrap::init_code;
using reflectomiter::io::iwrap::clean_up;
using reflectomiter::io::iwrap::code_step;
using reflectomiter::io::iwrap::get_code_state;
using reflectomiter::io::iwrap::restore_code_state;
using reflectomiter::io::iwrap::get_timestamp_cpp;
```

---

## Optional CMake target for the C++ adapter

The C++ adapter can also be built through an optional CMake target:

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

Configure and build with:

```bash
cmake -S . -B build_iwrap_adapter \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_IWRAP_ADAPTER=ON

cmake --build build_iwrap_adapter
```

---

## Installing iWrap locally

For local development, install iWrap from source:

```bash
cd ~
git clone https://github.com/iterorganization/iWrap.git
cd iWrap

python -m pip install --upgrade pip setuptools wheel setuptools_scm
python -m pip install -v -e .
```

Check the installation:

```bash
iwrap --list-actor-types
```

If `pip install iwrap` is not available, use the source installation above.

---

## IMAS Access Layer visibility

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

On ITER or HPC systems, this should normally be handled by the official IMAS/iWrap environment setup or modules.

---

## C++ MUSCLE3/libmuscle requirements

Generating a full compiled C++ iWrap/MUSCLE3 actor may additionally require:

```bash
pkg-config --modversion ymmsl
pkg-config --modversion libmuscle_mpi
```

and the C++ header:

```cpp
#include <libmuscle/libmuscle.hpp>
```

If these are not available on WSL, use the Python wrapper locally and build the full compiled actor on an ITER/Gateway environment where IMAS, iWrap, MUSCLE3, MPI, and compiler modules are provided.

---

## Generate the iWrap actor

Once iWrap is installed, generate the actor with:

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

The generated actor is installed by default under:

```text
$HOME/IWRAP_ACTORS/
```

The current generated actor path still launches the external `reflectomITER` executable. It is not yet a true linked solver-library actor.

---

## Recommended workflow

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
