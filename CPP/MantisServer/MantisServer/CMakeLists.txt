cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -lpthread")
set(CMAKE_BUILD_TYPE Debug)

project(MantisServer)

set(PROJECT_LIBS "")

#find_package(HighFive REQUIRED)
#find_package(HDF5 REQUIRED)

find_package (Threads REQUIRED)

find_package(Boost 1.69 REQUIRED COMPONENTS serialization program_options system)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ) #${HDF5_INCLUDE_DIRS}
set(PROJECT_LIBS ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ) #${HDF5_C_LIBRARIES}


add_executable(${PROJECT_NAME} "MantisServer.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${PROJECT_LIBS})