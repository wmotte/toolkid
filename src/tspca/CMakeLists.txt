project( tspca )

include_directories( ../cmdparser )

add_executable( tspca main.cpp tspca.h tspca.cpp )
target_link_libraries( tspca cmdparser ITKCommon ITKIO ITKNumerics ITKStatistics )
