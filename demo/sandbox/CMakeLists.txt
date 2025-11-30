cmake_minimum_required(VERSION 3.16)
project(eugene-sandbox)

set(SOURCES sbox/layer.cpp sbox/main.cpp)

add_executable(eugene-sandbox ${SOURCES})

target_link_libraries(eugene-sandbox PRIVATE eugene-sph onyx)
