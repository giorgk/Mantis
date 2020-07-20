cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -lpthread")
set(CMAKE_BUILD_TYPE Release)

project(MantisServer)

find_package( CGAL QUIET COMPONENTS  )
include( ${CGAL_USE_FILE} )

find_package (Threads REQUIRED)

find_package(Boost 1.69 REQUIRED COMPONENTS serialization program_options system)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})


add_executable(${PROJECT_NAME} "MantisServer.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${CGAL_LIBRARIES})