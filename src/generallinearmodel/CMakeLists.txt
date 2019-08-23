cmake_minimum_required(VERSION 2.6)
project( generallinearmodel )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../nr ../fftw++ )

if( Boost_FOUND )

  # General Linear model with autocorrelation prewhitening 
  add_executable( generallinearmodel generallinearmodel.cpp )
  target_link_libraries( generallinearmodel cmdparser fftwpp fftw3 ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
