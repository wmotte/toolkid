project( fdfconvert )

include_directories( ../cmdparser )

add_executable( fdfconvert fdfConvert.h fdfDictionary.h fdfConvert.cpp main.cpp )
target_link_libraries( fdfconvert cmdparser ITKCommon ITKIO )
