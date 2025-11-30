cmake_minimum_required(VERSION 3.16)
project(eugene)

set(SPH_SOURCES gene/sph/kernels.cpp)

add_library(eugene-sph STATIC ${SPH_SOURCES})

include(FetchContent)
FetchContent_Declare(
  toolkit
  GIT_REPOSITORY https://github.com/ismawno/toolkit.git
  GIT_TAG main
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(toolkit)

tkit_register_for_reflection(eugene-sph SOURCES gene/sph/kernels.hpp)
tkit_register_for_yaml_serialization(eugene-sph SOURCES gene/sph/kernels.hpp)

target_include_directories(eugene-sph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_precompile_headers(eugene-sph PRIVATE
                          ${CMAKE_CURRENT_SOURCE_DIR}/gene/core/pch.hpp)

target_link_libraries(eugene-sph PUBLIC toolkit)
target_compile_definitions(
  eugene-sph PUBLIC GENE_ROOT_PATH="${EUGENE_ROOT_PATH}"
                    GENE_VERSION=\"v0.1.0\")

tkit_default_configure(eugene-sph)
