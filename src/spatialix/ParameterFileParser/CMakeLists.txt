project( ParameterFileParser )

add_library(
 ParameterFileParser
 itkParameterFileParser.h
 itkParameterFileParser.cxx
 itkParameterMapInterface.h
 itkParameterMapInterface.cxx
)

target_link_libraries( ParameterFileParser ITKCommon )

