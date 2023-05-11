cmake_minimum_required(VERSION 3.15)

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -lpthread")
endif(UNIX)
##set(CMAKE_CXX_STANDARD 11)
##set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

project(MantisClient)

find_package (Threads REQUIRED)

find_package(Boost 1.69 REQUIRED serialization thread system)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})


add_executable(${PROJECT_NAME} "TestClient.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})