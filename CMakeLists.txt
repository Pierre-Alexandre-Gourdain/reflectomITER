cmake_minimum_required(VERSION 3.25)

option(USE_CUDA "Enable CUDA GPU support" ON)
option(USE_IMAS "Enable IMASH to use IMAS input" ON)
option(BUILD_DOCS "Build Doxygen documentation" ON)


if(USE_CUDA)
    project(reflectomITER LANGUAGES C CXX CUDA Fortran)
else()
    project(reflectomITER LANGUAGES C CXX Fortran)
endif()

# ------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------
set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------
# CUDA configuration
# ------------------------------------------------------------
if(USE_CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES 86 CACHE STRING "CUDA architectures")
endif()

# -------------------------------------------------------------
# Compile definitions
# -------------------------------------------------------------
add_compile_definitions(FD_ORDER=2)

# ------------------------------------------------------------
# AMReX configuration
# ------------------------------------------------------------

set(AMREX_USE_MPI ON CACHE BOOL "" FORCE)

if(USE_CUDA)
    set(AMREX_USE_CUDA ON CACHE BOOL "" FORCE)
else()
    set(AMREX_USE_CUDA OFF CACHE BOOL "" FORCE)
endif()

# First try to find AMReX normally:
#   - from CMAKE_PREFIX_PATH
#   - from AMReX_DIR if supplied on the command line
#   - from system/module paths
find_package(AMReX QUIET)

if(USE_CUDA)
	# If not found, fall back to the default local GPU installation
	if(NOT AMReX_FOUND)
		message(STATUS "AMReX not found in default CMake paths.")
		message(STATUS "Falling back to $ENV{HOME}/amrex_GPU/installdir/lib/cmake/AMReX")

		set(AMReX_DIR
			"$ENV{HOME}/amrex_GPU/installdir/lib/cmake/AMReX"
			CACHE PATH "Path to AMReX CMake config"
			FORCE
		)

		find_package(AMReX REQUIRED)
	endif()
else()
	# If not found, fall back to the default local CPU installation
	if(NOT AMReX_FOUND)
		message(STATUS "AMReX not found in default CMake paths.")
		message(STATUS "Falling back to $ENV{HOME}/amrex_no_GPU/installdir/lib/cmake/AMReX")

		set(AMReX_DIR
			"$ENV{HOME}/amrex_no_GPU/installdir/lib/cmake/AMReX"
			CACHE PATH "Path to AMReX CMake config"
			FORCE
		)

		find_package(AMReX REQUIRED)
	endif()
endif()
# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------
find_package(MPI REQUIRED)

# ------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------
# IMASH (optional)
# ------------------------------------------------------------
if(USE_IMAS)

    include(FetchContent)

    set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/_deps)

    set(IMASH_BUILD_EXECUTABLE OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        IMASH
        GIT_REPOSITORY https://github.com/Pierre-Alexandre-Gourdain/IMASH.git
        GIT_TAG main
    )

    FetchContent_MakeAvailable(IMASH)

    if(NOT TARGET IMASH)
        message(FATAL_ERROR "IMASH was fetched, but target IMASH was not created.")
    endif()

endif()

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
file(GLOB_RECURSE FDTD_SOURCES CONFIGURE_DEPENDS
    src/*.cpp
)

if(NOT USE_IMAS)
    list(FILTER FDTD_SOURCES EXCLUDE REGEX "ImasReader.cpp")
endif()

add_executable(reflectomITER ${FDTD_SOURCES})

target_include_directories(reflectomITER PRIVATE include)

# ------------------------------------------------------------
# CUDA source handling
# ------------------------------------------------------------
if(USE_CUDA AND AMReX_CUDA)
    set_target_properties(reflectomITER PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )

    target_compile_options(reflectomITER PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda --expt-relaxed-constexpr -Xcompiler=-fopenmp>
    )

    set_source_files_properties(
        src/fdtd/EHSolver.cpp
        src/fdtd/EHJSolver.cpp
        src/fdtd/PMLSolver.cpp
        src/fdtd/DivCleaning.cpp
        src/main.cpp
        src/io/Config.cpp
        src/utils/HelperFunctions.cpp
        src/sources/PointSource.cpp
        src/sources/ADEPlasma.cpp
        PROPERTIES LANGUAGE CUDA
    )
endif()

#-------------------------------------------------------------
# iWrap handling
#-------------------------------------------------------------
option(BUILD_IWRAP_ADAPTER "Build reflectomITER iWrap adapter" OFF)

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
# ------------------------------------------------------------
# Linking
# ------------------------------------------------------------
target_link_libraries(reflectomITER PRIVATE
    AMReX::amrex
    MPI::MPI_CXX
	OpenMP::OpenMP_CXX
)

if(USE_IMAS)
    target_link_libraries(reflectomITER PRIVATE IMASH)
endif()

if(USE_IMAS)
    target_compile_definitions(reflectomITER PRIVATE USE_IMAS=1)
endif()

# ------------------------------------------------------------
# Doxygen documentation
#     to build do 
#          cmake -S . -B build
#          cmake --build build --target docs
# ------------------------------------------------------------
option(BUILD_DOCS "Build Doxygen documentation" ON)

# ------------------------------------------------------------
# Doxygen documentation
# ------------------------------------------------------------
option(BUILD_DOCS "Build Doxygen documentation" ON)

if(BUILD_DOCS)

    find_package(Doxygen)

    if(DOXYGEN_FOUND)

        message(STATUS "Doxygen found: ${DOXYGEN_EXECUTABLE}")

        set(DOXYGEN_PROJECT_NAME "reflectomITER")
        set(DOXYGEN_PROJECT_BRIEF "AMReX-based FDTD Maxwell solver")
        set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/docs")

        set(DOXYGEN_GENERATE_HTML YES)
        set(DOXYGEN_GENERATE_LATEX NO)
        set(DOXYGEN_GENERATE_TREEVIEW YES)

        # Math rendering
        set(DOXYGEN_USE_MATHJAX YES)
        set(DOXYGEN_MATHJAX_RELPATH "https://cdn.jsdelivr.net/npm/mathjax@2")
        set(DOXYGEN_FORMULA_FONTSIZE 10)

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

        if(EXISTS "${CMAKE_SOURCE_DIR}/src/main_page.md")
            set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_SOURCE_DIR}/src/main_page.md")

            doxygen_add_docs(
                docs
                "${CMAKE_SOURCE_DIR}/src/main_page.md"
                "${CMAKE_SOURCE_DIR}/include"
                "${CMAKE_SOURCE_DIR}/src"
                COMMENT "Generating Doxygen documentation"
            )
        else()
            doxygen_add_docs(
                docs
                "${CMAKE_SOURCE_DIR}/include"
                "${CMAKE_SOURCE_DIR}/src"
                COMMENT "Generating Doxygen documentation"
            )
        endif()

    else()
        message(WARNING "Doxygen not found: documentation target 'docs' will not be created")
    endif()

endif()
# ------------------------------------------------------------
# Build summary
# ------------------------------------------------------------
message(STATUS "---------------------------------------------")
message(STATUS "Build type:          ${CMAKE_BUILD_TYPE}")
message(STATUS "USE_CUDA:            ${USE_CUDA}")
message(STATUS "AMReX CUDA enabled:  ${AMReX_CUDA}")
message(STATUS "USE_IMAS:            ${USE_IMAS}")
message(STATUS "BUILD_DOCS:          ${BUILD_DOCS}")
message(STATUS "C++ Compiler:        ${CMAKE_CXX_COMPILER}")
message(STATUS "CUDA compiler:       ${CMAKE_CUDA_COMPILER}")
message(STATUS "MPI compiler:        ${MPI_CXX_COMPILER}")
message(STATUS "---------------------------------------------")