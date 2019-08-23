cmake_minimum_required( VERSION 2.6 )
project( bayesian_dti )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser )

if( Boost_FOUND )

  # Bayesian DTI Gibbs sampler
  add_executable( bayesian_dti bayesianDTI.cpp )
  target_link_libraries( bayesian_dti cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO  )

endif( Boost_FOUND )
