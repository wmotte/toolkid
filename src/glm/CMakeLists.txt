cmake_minimum_required(VERSION 2.6)
project( generalized_linear_model )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser )

if( Boost_FOUND )

  # Generalized Linear model
  add_executable( glm glm.cpp )
  target_link_libraries( glm cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
