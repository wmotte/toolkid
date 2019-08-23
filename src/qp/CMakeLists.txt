cmake_minimum_required(VERSION 2.6)

project( qp )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser )

add_library( qp
    Array.hh
    QuadProg++.hh
    Array.cc
    QuadProg++.cc
)
  
# Test quadratic programming
add_executable( test_qp main.cpp )
target_link_libraries( test_qp qp )
