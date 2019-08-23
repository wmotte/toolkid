project( spatialix )

include_directories( ../cmdparser Transforms ParameterFileParser)

subdirs(
  ParameterFileParser
)

add_library(
  spatialistix
  spatialix.h
  spatialix.cpp
  tkdTransform.h
  tkdTransform.cpp
  tkdBSplineTransform.h
  tkdBSplineTransform.cpp
  tkdTransformFactory.h
  tkdTransformFactory.cpp
  tkdCombinationTransform.h
  tkdCombinationTransform.cpp
  tkdSimilarityTransform.h
  tkdSimilarityTransform.cpp
  tkdAffineTransform.h
  tkdAffineTransform.cpp
)  

add_executable(
  spatialix
  main.cpp
)

target_link_libraries(
  spatialix
  cmdparser
  spatialistix
  ParameterFileParser
  ITKCommon
  ITKIO
)
