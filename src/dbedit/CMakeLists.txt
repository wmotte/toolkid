project( dbedit )

include_directories( include src ${Boost_INCLUDE_DIRS} ../cmdparser ../zlib )

# dbedit 
add_executable( dbedit dbedit.cpp zfstream.cc zfstream.h anautils.c anautils.h dbh.h )
target_link_libraries( dbedit cmdparser zlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

