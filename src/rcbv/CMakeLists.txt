project( rcbv )

include_directories( ../cmdparser )

add_executable( rcbv rcbv.cpp )
target_link_libraries( rcbv cmdparser ITKCommon ITKIO ITKNumerics )

