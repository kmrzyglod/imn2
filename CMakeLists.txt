cmake_minimum_required(VERSION 3.3)
project(imn2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp PotentialMesh.h)
add_executable(imn2 ${SOURCE_FILES})