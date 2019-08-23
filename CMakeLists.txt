cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
project( toolkid )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )
#include( UseITK.cmake )

# cern ROOT requirements
set( USE_ROOT ON CACHE BOOL "Use ROOT support" )
if( USE_ROOT )
    set( ROOT_HOME "" CACHE FILEPATH "ROOT home" )
    include( FindROOT.cmake )
endif( USE_ROOT )

# GSL requirements
set( USE_GSL ON CACHE BOOL "Use GSL support" )
if( USE_GSL )
    set( GSL_HOME "" CACHE FILEPATH "GSL home" )
    include( FindGSL.cmake )
endif( USE_GSL )

# Boost
set( Boost_USE_STATIC_LIBS ON )
set( Boost_USE_MULTITHREAD ON )

find_package( Boost COMPONENTS program_options system graph regex REQUIRED )

message( "Boost program_options library: " ${Boost_PROGRAM_OPTIONS_LIBRARY} )
message( "Boost system library: " ${Boost_SYSTEM_LIBRARY} )
message( "Boost graph library: " ${Boost_GRAPH_LIBRARY} )
message( "Boost regex library: " ${Boost_REGEX_LIBRARY} )

message( "Boost include directories: " ${Boost_INCLUDE_DIRS} )

set( USE_BLAS ON CACHE BOOL "Build with BLAS support" )
mark_as_advanced( USE_BLAS )

set( LIBRARY_OUTPUT_PATH ${toolkid_BINARY_DIR}/lib CACHE INTERNAL "Single library output directory" )
set( EXECUTABLE_OUTPUT_PATH ${toolkid_BINARY_DIR}/bin CACHE INTERNAL "Single binary output directory" )

subdirs( src )
