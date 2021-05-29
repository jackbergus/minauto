cmake_minimum_required(VERSION 3.16)
project(complete_matching)

set(CMAKE_CXX_STANDARD 11)


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost REQUIRED COMPONENTS filesystem regex system thread date_time wave)
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories(include/minauto)
include_directories(submodules/mpreal)

add_library(equiv source/BLAS.cpp source/Util.cpp source/Automaton.cpp source/WordTree.cpp source/RandomValueMap.cpp source/Timer.cpp #source/mpreal.cpp
        source/dlmalloc.cpp)

add_executable(minauto source/main.cpp)
target_link_libraries(minauto equiv mpfr gmp ${Boost_LIBRARIES})