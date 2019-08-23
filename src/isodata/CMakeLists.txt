cmake_minimum_required(VERSION 2.6)
project( isodata )

include_directories( ../cmdparser )

include_directories( ../kmeans )

# fast isodata 
add_executable( isodata isodata.cpp 
  IsoDataImage.cpp
  Point.cpp
  Point.h
)

target_link_libraries( 
  isodata 
  kmeans 
  cmdparser 
  ITKCommon 
  ITKIO 
)
