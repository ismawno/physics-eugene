cmake_minimum_required(VERSION 3.16)
project(eugene)

set(SOURCES gene/sph/kernels.cpp)

add_library(eugene STATIC ${SOURCES})

include(FetchContent)
FetchContent_Declare(
  toolkit
  GIT_REPOSITORY https://github.com/ismawno/toolkit.git
  GIT_TAG main
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(toolkit)

tkit_register_for_reflection(eugene SOURCES gene/sph/kernels.hpp)
tkit_register_for_yaml_serialization(eugene SOURCES gene/sph/kernels.hpp)

target_include_directories(eugene PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(eugene PRIVATE toolkit)
target_compile_definitions(eugene PRIVATE GENE_ROOT_PATH="${EUGENE_ROOT_PATH}"
                                          GENE_VERSION=\"v0.4.0\")
tkit_default_configure(eugene)
