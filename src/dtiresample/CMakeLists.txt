cmake_minimum_required(VERSION 2.6)
project( dtiresample )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser )

#Set( dtiProcessFiles
#    dtiprocessFiles/deformationfieldio.h
#    dtiprocessFiles/deformationfieldio.cxx
#    dtiprocessFiles/itkHFieldToDeformationFieldImageFilter.h
#    dtiprocessFiles/itkHFieldToDeformationFieldImageFilter.txx
#)

#add_library( dtiprocess ${dtiProcessFiles} )


if( Boost_FOUND )

  # dtiresample
  #add_executable( dtiresample dtiresample.cpp )
  #target_link_libraries( dtiresample dtiprocess cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO ITKStatistics ITKBasicFilters )

endif( Boost_FOUND )
