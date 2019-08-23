project( orient )

include_directories( ../cmdparser )

add_executable(
  orient
  orient.cpp
)

target_link_libraries(
  orient
  cmdparser
  ITKCommon
  ITKIO
)
