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

The next step is to generate the `reflectomITER` input file automatically from actor parameters or IDS data.
