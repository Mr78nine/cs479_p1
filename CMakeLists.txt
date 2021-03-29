cmake_minimum_required(VERSION 3.16)
project(cs479_proj1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS  "-lboost_iostreams -lboost_system -lboost_filesystem")
add_executable(cs479_proj1 main.cpp)