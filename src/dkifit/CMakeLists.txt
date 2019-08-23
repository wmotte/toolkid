project( dkifit )

include_directories( ../cmdparser ${Boost_INCLUDE_DIRS} ../uqp )

add_library( dkifit dkifitCommon.cpp dkiOptimizeFunction.cpp )

add_executable( dkifit_linear dkifit_linear.cpp )
target_link_libraries( dkifit_linear cmdparser dkifit ITKCommon ITKIO ${Boost_LIBRARIES} )

add_executable( dkifit_constrained dkifit_constrained.cpp )
target_link_libraries( dkifit_constrained cmdparser dkifit ITKCommon ITKIO ${Boost_LIBRARIES} )

add_executable( dkifit_nonlinear dkifit_nonlinear.cpp )
target_link_libraries( dkifit_nonlinear cmdparser dkifit ITKCommon ITKIO ${Boost_LIBRARIES} )

add_executable( dkifit_RESTORE dkifit_RESTORE.cpp )
target_link_libraries( dkifit_RESTORE cmdparser dkifit ITKCommon ITKIO ${Boost_LIBRARIES} )



add_executable( dkimaps dkimaps.cpp )
target_link_libraries( dkimaps cmdparser ITKCommon ITKIO ${Boost_LIBRARIES} )


add_executable( dki_test_index TestIndex.cpp )
target_link_libraries( dki_test_index ITKCommon ITKIO ${Boost_LIBRARIES} )