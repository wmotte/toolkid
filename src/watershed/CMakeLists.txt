cmake_minimum_required(VERSION 2.6)
project( watershed )

include_directories( ../cmdparser )

add_executable( watershed watershed.cpp )
target_link_libraries( watershed cmdparser ITKCommon ITKIO ITKAlgorithms ITKBasicFilters)

