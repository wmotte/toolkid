project( bfloatnii )

include_directories(
  ../cmdparser
)

add_executable( bfloat2nii bfloat2nii.cpp )
target_link_libraries( bfloat2nii cmdparser ITKCommon ITKIO ITKBasicFilters )

add_executable( nii2bfloat nii2bfloat.cpp )
target_link_libraries( nii2bfloat cmdparser ITKCommon ITKIO )

