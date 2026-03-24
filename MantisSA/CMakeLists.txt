cmake_minimum_required(VERSION 3.16)

project(MantisSA LANGUAGES CXX)

# -----------------------------
# User options
# -----------------------------
option(USEHF "Compile with HDF5 support" OFF)
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
# On Windows with vcpkg/MSVC, imported config packages usually work well.
# On Linux, older Boost installs may only work through the classic module mode.
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_NO_WARN_NEW_VERSIONS ON)

# Optional: this often helps on Linux with older system Boost installs
# while still working on Windows/vcpkg.
set(Boost_NO_BOOST_CMAKE ON)

find_package(Boost REQUIRED COMPONENTS mpi serialization program_options system)

target_include_directories(MantisSA PRIVATE ${Boost_INCLUDE_DIRS})

# Prefer imported targets if they exist
if(TARGET Boost::mpi)
    target_link_libraries(MantisSA PRIVATE
            Boost::mpi
            Boost::serialization
            Boost::program_options
            Boost::system
    )
else()
    # Fallback for older FindBoost behavior
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

    find_package(HDF5 REQUIRED COMPONENTS CXX)

    # HighFive may not exist on older Linux systems.
    # Use it only if available, otherwise continue with plain HDF5.
    find_package(HighFive QUIET)

    if(HighFive_FOUND)
        message(STATUS "HighFive found")
        if(TARGET HighFive)
            target_link_libraries(MantisSA PRIVATE HighFive)
        elseif(TARGET HighFive::HighFive)
            target_link_libraries(MantisSA PRIVATE HighFive::HighFive)
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
