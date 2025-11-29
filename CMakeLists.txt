cmake_minimum_required(VERSION 3.16...3.31)
project(eugene-root)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE
    "Debug"
    CACHE STRING "")

set(EUGENE_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(eugene)
