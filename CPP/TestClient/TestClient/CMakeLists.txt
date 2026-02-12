cmake_minimum_required(VERSION 3.15)

project(MantisClient LANGUAGES CXX)

# Portable C++ standard selection (do NOT hardcode -std=... in flags)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package (Threads REQUIRED)

# Boost with components
find_package(Boost REQUIRED COMPONENTS serialization thread system)

add_executable(${PROJECT_NAME} TestClient.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE
        Threads::Threads
        Boost::serialization
        Boost::thread
        Boost::system
)