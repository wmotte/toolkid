project( dkisim )

# Simulate crossing-fiber tensors for DKI validation study.

include_directories( ../cmdparser ../dkifit ${Boost_INCLUDE_DIRS} )

add_executable( dkisim dkisim.cpp )
target_link_libraries( dkisim cmdparser dkifit ITKCommon ITKIO ${Boost_LIBRARIES} )
