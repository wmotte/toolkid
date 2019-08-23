project( meants )

include_directories( ../cmdparser )

add_executable( 
  meants
  meants.cpp
)
  
target_link_libraries(
  meants
  cmdparser
  ITKCommon
  ITKIO
)
