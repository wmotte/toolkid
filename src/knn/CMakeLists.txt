project( knn )

include_directories( include src include isodata_src ${Boost_INCLUDE_DIRS} ../cmdparser ../flens ../graphs )

# Sources of non-templated classes.
  set(ISO_SRCS
      isodata_src/Point.h
      isodata_src/Point.cpp 
      isodata_src/IsoData.h
      isodata_src/IsoData.cpp
      isodata_src/KCtree.h
      isodata_src/KCtree.cpp
      isodata_src/KCutil.h
      isodata_src/KCutil.cpp
      isodata_src/KM_ANN.h
      isodata_src/KM_ANN.cpp
      isodata_src/KMcenters.h
      isodata_src/KMcenters.cpp
      isodata_src/KMdata.h
      isodata_src/KMdata.cpp
      isodata_src/KMrand.h
      isodata_src/KMrand.cpp
      isodata_src/KMeans.h
      isodata_src/KMeans.cpp
      isodata_src/KMfilterCenters.h  
      isodata_src/KMfilterCenters.cpp
  )

# create the ISOlib library.
add_library( ISOlib ${ISO_SRCS})

# Sources of non-templated classes.
  set(ANN_SRCS
      src/ANN.cpp
      src/bd_fix_rad_search.cpp
      src/bd_pr_search.cpp
      src/bd_search.cpp
      src/bd_tree.cpp
      src/brute.cpp
      src/kd_dump.cpp
      src/kd_fix_rad_search.cpp
      src/kd_pr_search.cpp
      src/kd_search.cpp
      src/kd_split.cpp
      src/kd_tree.cpp
      src/kd_util.cpp
      src/perf.cpp
     )

# create the ANNlib library.
add_library( ANNlib ${ANN_SRCS})

# sources of comman ann library.
  set( ANN_COMMON_SRCS
      annCommon.h
      annCommon.cpp
      annRelabel.h
      annRelabel.cpp
      itkLabelGeometryImageFilter.h
      itkLabelGeometryImageFilter.txx
     )

# Create ann common library
add_library( ann ${ANN_COMMON_SRCS} )

# bootstrap
add_executable( bootstrap bootstrap.cpp )
target_link_libraries( bootstrap cmdparser graph ann ANNlib ${Boost_LIBRARIES} ITKCommon ITKIO ITKStatistics ITKBasicFilters )



# clustering
add_executable( ann_clustering annClustering.cpp )
target_link_libraries( ann_clustering ann ANNlib graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO ITKStatistics )

# find optimal clustering
add_executable( ann_find_optimal_clustering annFindOptimalClustering.cpp )
target_link_libraries( ann_find_optimal_clustering ann ANNlib graph cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO ITKStatistics )

# modularity
add_executable( ann_modularity annModularity.cpp )
target_link_libraries( ann_modularity ann ANNlib cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO ITKStatistics )

# relabel
add_executable( ann_relabel relabel.cpp )
target_link_libraries( ann_relabel cmdparser ann ${Boost_LIBRARIES} ITKCommon ITKIO ITKStatistics )

# max
add_executable( ann_max annMax.cpp )
target_link_libraries( ann_max cmdparser ${Boost_LIBRARIES} )

# isodata 
add_executable( isodata isodata.cpp )
target_link_libraries( isodata ann cmdparser ANNlib ISOlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# ann border distance
  add_executable( ann_border_distance annBorderDistance.cpp )
target_link_libraries( ann_border_distance ann cmdparser ANNlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# ann2fig
  add_executable( ann2fig ann2fig.cpp )
target_link_libraries( ann2fig cmdparser ANNlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# ann train
  add_executable( ann_train annTrain.cpp )
target_link_libraries( ann_train ann ANNlib cmdparser ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# ann classify
  add_executable( ann_classify annClassify.cpp )
target_link_libraries( ann_classify ann cmdparser ANNlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# ann similarity
  add_executable( ann_similarity annSimilarity.cpp )
target_link_libraries( ann_similarity ann cmdparser ANNlib ITKCommon ITKIO ITKBasicFilters ITKStatistics ${Boost_LIBRARIES} )

# ann stats
  add_executable( ann_stats annStats.cpp )
target_link_libraries( ann_stats ann cmdparser ANNlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# ann minimal connected component
  add_executable( ann_min_connected_component annMinConnectedComponent.cpp )
target_link_libraries( ann_min_connected_component ann cmdparser ANNlib ITKCommon ITKIO ITKStatistics ${Boost_LIBRARIES} )

# example ann
#add_executable( ann_example ann.cpp )
#target_link_libraries( ann_example cmdparser ANNlib ITKCommon ITKIO ITKStatistics )

