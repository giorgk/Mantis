cmake_minimum_required(VERSION 3.16)

project(MantisSA LANGUAGES CXX)

# -----------------------------
# User options
# -----------------------------
option(USEHF "Compile with HDF5 support" ON)
option(USEZLIB "Compile with ZLIB compression support" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------
# Compiler warnings
# -----------------------------
if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# -----------------------------
# Executable
# -----------------------------
add_executable(MantisSA main.cpp)

# -----------------------------
# MPI
# -----------------------------
find_package(MPI REQUIRED COMPONENTS CXX)

if(TARGET MPI::MPI_CXX)
    target_link_libraries(MantisSA PRIVATE MPI::MPI_CXX)
else()
    # Fallback for older CMake / MPI packages
    target_include_directories(MantisSA PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(MantisSA PRIVATE ${MPI_CXX_LIBRARIES})
    if(MPI_CXX_COMPILE_OPTIONS)
        target_compile_options(MantisSA PRIVATE ${MPI_CXX_COMPILE_OPTIONS})
    endif()
    if(MPI_CXX_LINK_FLAGS)
        string(REPLACE " " ";" MPI_CXX_LINK_FLAGS_LIST "${MPI_CXX_LINK_FLAGS}")
        target_link_options(MantisSA PRIVATE ${MPI_CXX_LINK_FLAGS_LIST})
    endif()
endif()

# -----------------------------
# Boost
# -----------------------------
# On Windows with vcpkg/MSVC, config packages usually work well.
# On Linux, older Boost installs may still prefer classic module mode.
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_NO_WARN_NEW_VERSIONS ON)

# Keep old behavior compatible with non-vcpkg Linux installs
set(Boost_NO_BOOST_CMAKE ON)

find_package(Boost REQUIRED COMPONENTS mpi serialization program_options system)

if(TARGET Boost::mpi)
    target_link_libraries(MantisSA PRIVATE
            Boost::mpi
            Boost::serialization
            Boost::program_options
            Boost::system
    )
else()
    target_include_directories(MantisSA PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(MantisSA PRIVATE ${Boost_LIBRARIES})
endif()

# -----------------------------
# ZLIB (optional)
# -----------------------------
if(USEZLIB)
    message(STATUS "Compile with ZLIB compression")
    find_package(ZLIB REQUIRED)
    target_compile_definitions(MantisSA PRIVATE _USEZLIB=1)

    if(TARGET ZLIB::ZLIB)
        target_link_libraries(MantisSA PRIVATE ZLIB::ZLIB)
    else()
        target_include_directories(MantisSA PRIVATE ${ZLIB_INCLUDE_DIRS})
        target_link_libraries(MantisSA PRIVATE ${ZLIB_LIBRARIES})
    endif()
else()
    message(STATUS "Compile without ZLIB compression")
    target_compile_definitions(MantisSA PRIVATE _USEZLIB=0)
endif()

# -----------------------------
# HDF5 / HighFive (optional)
# -----------------------------
if(USEHF)
    message(STATUS "Compile with HDF5")
    target_compile_definitions(MantisSA PRIVATE _USEHF=1)

    # First try config-mode package discovery (works with vcpkg)
    find_package(hdf5 CONFIG QUIET)

    if(hdf5_FOUND)
        message(STATUS "Using HDF5 from package config")

        # HighFive from vcpkg/config mode if available
        find_package(HighFive CONFIG QUIET)
        if(NOT HighFive_FOUND)
            find_package(HighFive QUIET)
        endif()

        if(HighFive_FOUND)
            message(STATUS "HighFive found")
            if(TARGET HighFive::HighFive)
                target_link_libraries(MantisSA PRIVATE HighFive::HighFive)
            elseif(TARGET HighFive_HighFive)
                target_link_libraries(MantisSA PRIVATE HighFive_HighFive)
            elseif(TARGET HighFive)
                target_link_libraries(MantisSA PRIVATE HighFive)
            endif()
        else()
            message(STATUS "HighFive not found; continuing with HDF5 only")
        endif()

        # Link whichever HDF5 target exists
        if(TARGET hdf5::hdf5_cpp-shared)
            target_link_libraries(MantisSA PRIVATE hdf5::hdf5_cpp-shared)
        elseif(TARGET hdf5::hdf5_cpp)
            target_link_libraries(MantisSA PRIVATE hdf5::hdf5_cpp)
        elseif(TARGET hdf5::hdf5-shared)
            target_link_libraries(MantisSA PRIVATE hdf5::hdf5-shared)
        elseif(TARGET hdf5::hdf5)
            target_link_libraries(MantisSA PRIVATE hdf5::hdf5)
        else()
            message(FATAL_ERROR "HDF5 config package found, but no known hdf5 target is available")
        endif()

        # Optional HL library if present
        if(TARGET hdf5::hdf5_hl-shared)
            target_link_libraries(MantisSA PRIVATE hdf5::hdf5_hl-shared)
        elseif(TARGET hdf5::hdf5_hl)
            target_link_libraries(MantisSA PRIVATE hdf5::hdf5_hl)
        endif()

    else()
        message(STATUS "Falling back to CMake FindHDF5 module")
        find_package(HDF5 REQUIRED COMPONENTS CXX)

        # HighFive may not exist on older Linux systems.
        # Use it only if available, otherwise continue with plain HDF5.
        find_package(HighFive QUIET)

        if(HighFive_FOUND)
            message(STATUS "HighFive found")
            if(TARGET HighFive::HighFive)
                target_link_libraries(MantisSA PRIVATE HighFive::HighFive)
            elseif(TARGET HighFive)
                target_link_libraries(MantisSA PRIVATE HighFive)
            endif()
        else()
            message(STATUS "HighFive not found; continuing with HDF5 only")
        endif()

        if(TARGET HDF5::HDF5)
            target_link_libraries(MantisSA PRIVATE HDF5::HDF5)
        else()
            target_include_directories(MantisSA PRIVATE ${HDF5_INCLUDE_DIRS})
            target_link_libraries(MantisSA PRIVATE ${HDF5_LIBRARIES})
        endif()
    endif()
else()
    message(STATUS "Compile without HDF5")
    target_compile_definitions(MantisSA PRIVATE _USEHF=0)
endif()

# -----------------------------
# Platform-specific notes
# -----------------------------
if(WIN32)
    message(STATUS "Configuring for Windows")
    target_compile_definitions(MantisSA PRIVATE NOMINMAX)
elseif(UNIX)
    message(STATUS "Configuring for Unix/Linux")
endif()