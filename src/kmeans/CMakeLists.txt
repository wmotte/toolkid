project( kmeans )

include_directories( ../cmdparser )

# sources of kmeans library.
set( KMEANS_SRCS
  src/DenseVector.cpp
  src/DenseVector.h
  src/Kmeans.cpp
  src/Kmeans.h
  src/main.cpp
  src/SimpleKmeans.cpp
  src/SimpleKmeans.h
  src/SparseVector.cpp
  src/SparseVector.h
  src/Utils.cpp
  src/Utils.h
  src/VectorFactory.cpp
  src/VectorFactory.h
  src/Vector.h
)

# Create kmeans common library
add_library( kmeans ${KMEANS_SRCS} )

add_executable( kmeanspp ${KMEANS_SRCS} )
target_link_libraries( kmeanspp kmeans cmdparser ITKCommon ITKIO )
