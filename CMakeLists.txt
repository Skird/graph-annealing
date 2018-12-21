cmake_minimum_required(VERSION 3.12)
project(graph_annealing)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -lboost_program_options)

add_executable(graph_annealing main.cc)
