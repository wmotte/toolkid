cmake_minimum_required(VERSION 2.6)
project( ap )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../graphs ../ann ../flens )

if( Boost_FOUND )

  # Affinity Propagation
  add_executable( affinity_propagation affinityPropagation.cpp )
  target_link_libraries( affinity_propagation cmdparser graph ann ITKCommon ITKIO ITKStatistics )

endif( Boost_FOUND )
