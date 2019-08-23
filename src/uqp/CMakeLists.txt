cmake_minimum_required(VERSION 2.6)

project( uqp )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser . )
  
# Test vnl and uBLAS quadratic programming
add_executable( test_uqp main.cpp )
target_link_libraries( test_uqp ITKCommon )
