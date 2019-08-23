cmake_minimum_required(VERSION 2.6)
project( list2mask )


include_directories( ../cmdparser )

add_executable( list2mask list2mask.cpp )
target_link_libraries( list2mask cmdparser ITKCommon ITKIO )
