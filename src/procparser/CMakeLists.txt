project( procparser )

include_directories( ../cmdparser )

add_library(
  procparser
  procparser.h
  procparser.cpp
)

add_executable( 
  pp 
  main.cpp
)

target_link_libraries(
  pp
  procparser
  cmdparser
)
