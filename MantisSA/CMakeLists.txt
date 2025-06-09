cmake_minimum_required(VERSION 3.27)
project(MantisSA)



set(CMAKE_CXX_STANDARD 17)

# MPI
find_package(MPI REQUIRED)

# Optional HDF5
if(USEHF)
    message(STATUS "Compile with HDF5")
    find_package(HighFive REQUIRED)
    find_package(HDF5 REQUIRED)
    add_definitions(-D_USEHF=1)
else()
    message(STATUS "Compile without HDF5")
    add_definitions(-D_USEHF=0)
endif(USEHF)

# Optional ZLIB
if(USEZLIB)
    message(STATUS "Compile with ZLIB compression")
    find_package(ZLIB REQUIRED)
    add_definitions(-D_USEZLIB=1)
else()
    message(STATUS "Compile without ZLIB compression")
    add_definitions(-D_USEZLIB=0)
endif()

# Boost
find_package(Boost 1.58 REQUIRED COMPONENTS mpi serialization program_options system)

# Include directories
include_directories(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
if(USEHF)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

# Executable
add_executable(MantisSA main.cpp)

# Link libraries
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${MPI_LIBRARIES})
if(USEHF)
    target_link_libraries(${PROJECT_NAME} ${HDF5_C_LIBRARIES})
endif()
if(USEZLIB)
    target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
endif()

# MPI compile/link flags
if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()