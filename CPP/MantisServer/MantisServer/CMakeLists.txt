cmake_minimum_required(VERSION 3.15)

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -lpthread")
endif(UNIX)

set(CMAKE_BUILD_TYPE Release)

project(MantisServer)

set(PROJECT_LIBS "")

if(USEHF)
    message(STATUS "Compile with HDF5")
    find_package(HighFive REQUIRED)
    find_package(HDF5 REQUIRED)
    add_definitions(-D_USEHF=1)
else()
    message(STATUS "Compile without HDF5")
    add_definitions(-D_USEHF=0)
endif(USEHF)

find_package (Threads REQUIRED)

find_package(Boost 1.69 REQUIRED COMPONENTS serialization program_options system)

if(USEHF)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
    set(PROJECT_LIBS ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${HDF5_C_LIBRARIES})
else()
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ) #${HDF5_INCLUDE_DIRS}
    set(PROJECT_LIBS ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ) #${HDF5_C_LIBRARIES}
endif(USEHF)


add_executable(${PROJECT_NAME} "MantisServer.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${PROJECT_LIBS})