cmake_minimum_required(VERSION 2.6)
project( em )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../flens )

if( Boost_FOUND )

  # expectation maximization
  add_executable( em em.cpp ) 
  target_link_libraries( em cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )

  # test123
  add_executable( test123 test.cpp ) 
  target_link_libraries( test123 cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO  )


endif( Boost_FOUND )
