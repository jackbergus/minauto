cmake_minimum_required(VERSION 3.16)
project(complete_matching)

set(CMAKE_CXX_STANDARD 20)



include_directories(include/minauto)
include_directories(submodules/mpreal)
include_directories(submodules/ublas/include)
include_directories(submodules/boost_core/include)

add_library(equiv source/BLAS.cpp source/Util.cpp source/Automaton.cpp source/WordTree.cpp source/RandomValueMap.cpp source/Timer.cpp #source/mpreal.cpp
        source/dlmalloc.cpp)

add_executable(minauto source/main.cpp)
target_link_libraries(minauto equiv mpfr gmp #${Boost_LIBRARIES}
        )