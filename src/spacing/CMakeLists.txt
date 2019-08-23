cmake_minimum_required(VERSION 2.6)
project( spacing )


include_directories( ../cmdparser )

# change spacing 
add_executable( changespacing changespacing.cpp )
target_link_libraries( changespacing cmdparser ITKCommon ITKIO )

# divide spacing
add_executable( dividespacing dividespacing.cpp )
target_link_libraries( dividespacing cmdparser ITKCommon ITKIO )

# set origin
add_executable( setorigin setorigin.cpp )
target_link_libraries( setorigin cmdparser ITKCommon ITKIO )

# get cog
add_executable( cog cog.cpp )
target_link_libraries( cog cmdparser ITKCommon ITKIO )
