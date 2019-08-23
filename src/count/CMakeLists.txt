cmake_minimum_required(VERSION 2.6)
project( voxelcount )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../flens )

if( Boost_FOUND )

  # count voxels above threshold
  add_executable( voxelcount voxelcount.cpp )
  target_link_libraries( voxelcount cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )
  
  add_executable( printcoordinates printcoordinates.cpp )
  target_link_libraries( printcoordinates cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
