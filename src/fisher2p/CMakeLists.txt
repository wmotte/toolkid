project( fisher2p )

include_directories( ../cmdparser )

add_executable( fisher2p fisher2p.cpp )
target_link_libraries( fisher2p cmdparser ITKCommon ITKIO ITKNumerics )

