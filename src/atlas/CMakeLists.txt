cmake_minimum_required(VERSION 2.6)
project( atlas )

include_directories( 
	../cmdparser 
	../procfid
  ../pp
  ../procparser
)

add_executable( sva2volume sva2volume.cpp )
target_link_libraries( sva2volume procparser cmdparser FID ITKCommon ITKIO ITKAlgorithms ITKBasicFilters)

add_executable( png2nii png2nii.cpp )
target_link_libraries( png2nii procparser cmdparser FID ITKCommon ITKIO ITKAlgorithms ITKBasicFilters)

