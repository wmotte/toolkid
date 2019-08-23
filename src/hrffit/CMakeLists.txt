cmake_minimum_required(VERSION 2.6)
project( hrffit )

include_directories( ../cmdparser )

add_executable(
  hrffit
  tkdHRFFit.h
  tkdHRFFit.cpp
  tkdGammaHRF.h
  tkdGammaHRF.cpp
  tkdFitGammaHRF.h
  tkdFitGammaHRF.cpp
  nrGamma.h
  nrGamma.cpp
  main.cpp
)

target_link_libraries(
  hrffit
  cmdparser
  ITKCommon
  ITKIO
)
