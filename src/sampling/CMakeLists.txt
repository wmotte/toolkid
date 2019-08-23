project( downsample )

include_directories( ../cmdparser )

add_executable( downsample downsample.cpp )
target_link_libraries( downsample cmdparser ITKCommon ITKIO ITKBasicFilters )

add_executable( downsample3d downsample3d.cpp )
target_link_libraries( downsample3d cmdparser ITKCommon ITKIO ITKBasicFilters )

add_executable( sample sample.cpp )
target_link_libraries( sample cmdparser ITKCommon ITKIO ITKBasicFilters )
