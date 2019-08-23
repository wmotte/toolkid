project( dtifit )

include_directories(
  ../cmdparser
  ../procparser
)

add_library(
  old.dtifit
  dtifit.h
  dtifit.cpp
)

add_executable(
  fitdti
  main.cpp
)

target_link_libraries(
  fitdti
  old.dtifit
  procparser
  cmdparser
  ITKCommon
  ITKIO
)
