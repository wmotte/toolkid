project( expandpng )

include_directories( ../cmdparser )

add_executable( expandpng expandpng.cpp )
target_link_libraries( expandpng cmdparser ITKCommon ITKIO )
