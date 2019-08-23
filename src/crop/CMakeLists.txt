cmake_minimum_required(VERSION 2.6)
project( crop )

include_directories( ../cmdparser )

# crop...
add_executable( crop crop.cpp )
target_link_libraries( crop cmdparser ITKCommon ITKIO )

add_executable( multicrop multicrop.cpp )
target_link_libraries( multicrop cmdparser ITKCommon ITKIO )

add_executable( cropf3 cropf3.cpp )
target_link_libraries( cropf3 cmdparser ITKCommon ITKIO )

add_executable( cropf3mask cropf3mask.cpp )
target_link_libraries( cropf3mask cmdparser ITKCommon ITKIO )

add_executable( clearregion clearregion.cpp )
target_link_libraries( clearregion cmdparser ITKCommon ITKIO )
