cmake_minimum_required(VERSION 2.6)

project( cqp )

include_directories( ${Boost_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ../cmdparser . )

add_library( CQP
    gsl_cqp.c
    initial_point.c
    mg_pdip.c
    ConstrainedQuadraticProgramming.cpp
)

add_executable( test_cqp main.cpp )

target_link_libraries( test_cqp CQP ${GSL_LIBRARIES} ITKCommon )
