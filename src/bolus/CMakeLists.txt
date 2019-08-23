cmake_minimum_required(VERSION 2.6)
project( bolus )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../flens )

if( Boost_FOUND )

  # bolus SVD
  add_executable( bolus_svd bolus_svd.cpp )
  target_link_libraries( bolus_svd cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
