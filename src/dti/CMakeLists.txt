project( dtifit )

include_directories(
  ../cmdparser
  ../procparser
)

add_library(
  dtifit
  dtifit.h
  dtifit.cpp
)

add_executable(
  dti
  main.cpp
)

target_link_libraries(
  dti
  dtifit
  procparser
  cmdparser
  ITKCommon
  ITKIO
)
