project( graphs )

include_directories( ${Boost_INCLUDE_DIRS} ../cmdparser ../flens ../nr )

if( Boost_FOUND )
	
	# graph library
  	add_library( graph 
  		graphCommon.h 
  		graphCommon.cpp 
  		Kmeans.h
  		Kmeans.cpp
  		Vector.h
  		DenseVector.h
  		DenseVector.cpp
  		SparseVector.h
  		SparseVector.cpp
  		Utils.h
  		Utils.cpp
  		VectorFactory.h
  		VectorFactory.cpp
  	)
	
	# test
	add_executable( graph_test test.cpp )
	target_link_libraries( graph_test flens blas graph cmdparser ${Boost_LIBRARIES} ITKIO ITKCommon )
	
	# discrete entropy
	add_executable( graph_entropy discreteEntropy.cpp )
	target_link_libraries( graph_entropy flens blas graph cmdparser ${Boost_LIBRARIES} ITKIO ITKCommon )
	
  	# laplacian
	add_executable( graph_laplacian laplacian.cpp )
	target_link_libraries( graph_laplacian flens blas graph cmdparser ${Boost_LIBRARIES} ITKIO ITKCommon )

	# segmentation combiner
	add_executable( graph_segmentation_combiner segmentationCombiner.cpp )
	target_link_libraries( graph_segmentation_combiner flens blas nr graph cmdparser ${Boost_LIBRARIES} ITKIO ITKCommon )

	# spectral matrix
	add_executable( graph_spectral_matrix spectralMatrix.cpp spectralMatrix.h )
	target_link_libraries( graph_spectral_matrix flens blas nr graph cmdparser ${Boos_LIBRARIES} ITKCommon ITKIO )

 # clustering coefficient
  add_executable( graph_clusteringcoefficient clusteringcoefficient.cpp )
  target_link_libraries( graph_clusteringcoefficient flens blas cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # threshold
  add_executable( graph_threshold threshold.cpp )
  target_link_libraries( graph_threshold graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # minimum spanning tree
  add_executable( graph_minimumspanningtree minimumspanningtree.cpp )
  target_link_libraries( graph_minimumspanningtree graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # synchronizability
  add_executable( graph_synchronizability synchronizability.cpp )
  target_link_libraries( graph_synchronizability cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # hierarchy
  add_executable( graph_hierarchy hierarchy.cpp )
  target_link_libraries( graph_hierarchy graph flens blas cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

	# list -> 3D mask
  	add_executable( graph_list2image list2image.cpp )
  	target_link_libraries( graph_list2image cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO )

	# betweenness centrality
  	add_executable( graph_betweenness betweenness.cpp )
  	target_link_libraries( graph_betweenness cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

	# create regular graph
	add_executable( graph_create_regular createRegular.cpp )
	target_link_libraries( graph_create_regular cmdparser ITKCommon ITKIO )

	# create smallworld graph
	add_executable( graph_create_smallworld createSmallWorld.cpp )
	target_link_libraries( graph_create_smallworld graph ${Boost_LIBRARIES} cmdparser ITKCommon ITKIO )


  # export graph
  add_executable( graph_export export.cpp )
  target_link_libraries( graph_export cmdparser graph ${Boost_LIBRARIES} ITKCommon ITKIO )

  # visability graph map
  add_executable( graph_visabilitygraphmap vgmap.cpp )
  target_link_libraries( graph_visabilitygraphmap graph ${Boost_LIBRARIES} ITKCommon ITKIO )

  # power law library
  add_library( powerlaw powerlawCommon.cpp )

  # power law fit
  add_executable( graph_powerlaw_fit powerlawFit.cpp )
  target_link_libraries( graph_powerlaw_fit powerlaw cmdparser ${Boost_LIBRARIES} )

  # degree
  add_executable( graph_degree degree.cpp )
  target_link_libraries( graph_degree graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # visability graph
  add_executable( graph_visabilitygraph visabilitygraph.cpp )
  target_link_libraries( graph_visabilitygraph graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # histogram
  add_executable( graph_histogram histogram.cpp )
  target_link_libraries( graph_histogram graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO )

  # imagehistogram
  add_executable( graph_imagehistogram imagehistogram.cpp )
  target_link_libraries( graph_imagehistogram ${Boost_LIBRARIES} cmdparser ITKCommon ITKIO ITKAlgorithms )

endif( Boost_FOUND )

