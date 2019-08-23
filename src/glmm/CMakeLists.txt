cmake_minimum_required(VERSION 2.6)
project( glmm )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../flens )

if( Boost_FOUND )

  # select GLMM voxels
  add_executable( glmm_select glmmSelectVoxels )
  target_link_libraries( glmm_select cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
