cmake_minimum_required(VERSION 2.6)
project( dynamics )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../flens )

if( Boost_FOUND )

  # logan (reversible model)
  add_executable( dynamics dynamics.cpp ) 
  target_link_libraries( dynamics cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
