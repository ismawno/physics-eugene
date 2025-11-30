cmake_minimum_required(VERSION 3.16)
project(eugene-demo)

include(FetchContent)
FetchContent_Declare(
  onyx
  GIT_REPOSITORY https://github.com/ismawno/onyx.git
  GIT_TAG main
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)

FetchContent_MakeAvailable(onyx)

add_subdirectory(sandbox)
