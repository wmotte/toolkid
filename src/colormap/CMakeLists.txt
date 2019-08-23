cmake_minimum_required(VERSION 2.6)
project( colormap )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../flens )

if( Boost_FOUND )

  # Colormap
  add_executable( colormap colormap.cpp )
  target_link_libraries( colormap cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
