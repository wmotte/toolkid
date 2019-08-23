project( mrfit )

include_directories( ../cmdparser ../procparser ../procfid )

add_library(
mrfit
itkMRT1ParameterMap3DImageFilter.h
itkMRT1ParameterMap3DImageFilter.txx
itkMRT2ParameterMap3DImageFilter.h
itkMRT2ParameterMap3DImageFilter.txx
mrfitCommon.h
mrfitCommon.cpp
looklocker.h
looklocker.cpp
)

add_executable(
t1map_looklocker
t1map_looklocker.cpp
)

add_executable(
t1map
t1map.cpp
)

add_executable(
t2map
t2map.cpp
)

target_link_libraries(
t1map
mrfit
cmdparser
procparser
ITKCommon
ITKIO
)

target_link_libraries(
t1map_looklocker
mrfit
FID
cmdparser
procparser
ITKCommon
ITKIO
)


target_link_libraries(
t2map
mrfit
cmdparser
procparser
ITKCommon
ITKIO
)

