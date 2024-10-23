cmake_minimum_required(VERSION 3.27)
project(MantisSA)



set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)

if(USEHF)
    message(STATUS "Compile with HDF5")
    find_package(HighFive REQUIRED)
    find_package(HDF5 REQUIRED)
    add_definitions(-D_USEHF=1)
else()
    message(STATUS "Compile without HDF5")
    add_definitions(-D_USEHF=0)
endif(USEHF)

find_package(Boost 1.58 REQUIRED COMPONENTS mpi serialization program_options system)
if(USEHF)
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
else()
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
endif(USEHF)

include_directories(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})

add_executable(MantisSA main.cpp)

if(USEHF)
    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${HDF5_C_LIBRARIES})
else()
    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${MPI_LIBRARIES})
endif(USEHF)

if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()