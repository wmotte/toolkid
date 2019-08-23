project( bandpass )

include_directories( ../cmdparser )

add_executable( bandpass
  bandpass.cpp
)

target_link_libraries( bandpass
  cmdparser
  ITKCommon
  ITKIO
  ITKNumerics
)

add_executable( spectrum
  spectrum.cpp
)

target_link_libraries( spectrum
  cmdparser
  ITKCommon
  ITKIO
  ITKNumerics
)

add_executable( removemean
  removemean.cpp
)

target_link_libraries( removemean
  cmdparser
  ITKCommon
  ITKIO
  ITKNumerics
)

add_executable( meanregression
  meanregression.cpp
)

target_link_libraries( meanregression
  cmdparser
  ITKCommon
  ITKIO
  ITKNumerics
)

add_executable( rsfilter
  rsfilter.cpp
)

target_link_libraries( rsfilter
  cmdparser
  ITKCommon
  ITKIO
  ITKNumerics
)

add_executable( linearregression
  linearregression.cpp
)

target_link_libraries( linearregression
  cmdparser
  ITKCommon
  ITKIO
  ITKNumerics
)
