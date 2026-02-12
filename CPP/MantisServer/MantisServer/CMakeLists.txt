cmake_minimum_required(VERSION 3.15)

project(MantisServer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package (Threads REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization program_options system)

add_executable(${PROJECT_NAME} "MantisServer.cpp")

target_link_libraries(${PROJECT_NAME} PRIVATE
        Threads::Threads
        Boost::serialization
        Boost::program_options
        Boost::system
)

if(USEHF)
    message(STATUS "Compile with HDF5")
    find_package(HighFive REQUIRED)
    find_package(HDF5 REQUIRED)
    add_definitions(-D_USEHF=1)
else()
    message(STATUS "Compile without HDF5")
    add_definitions(-D_USEHF=0)
endif(USEHF)

if(USEHF)
    target_link_libraries(${PROJECT_NAME} PRIVATE
            HighFive
            ${HDF5_C_LIBRARIES}
    )
    target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

