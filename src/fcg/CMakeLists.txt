project( fcg )

include_directories( ../cmdparser ../flens )

add_executable( rmatrix rmatrix.h rmatrix.cpp )
target_link_libraries( rmatrix cmdparser ITKCommon ITKIO )

add_library( rs
  rsCommon.h
  rsCommon.cpp
  rsConnectivityMatrix.h
  rsConnectivityMatrix.cpp
  rsPCA.h
  rsPCA.cpp
)

add_executable( indices
  indices.cpp
)

target_link_libraries( indices
  cmdparser
  ITKCommon
  ITKIO
)

add_executable( connectivitymatrix
  connectivitymatrix.cpp
)

target_link_libraries( connectivitymatrix
  rs
  cmdparser
  ITKCommon
  ITKIO
)

add_executable( matrix2csv
  matrix2csv.cpp
)

target_link_libraries( matrix2csv
  cmdparser
  ITKCommon
  ITKIO
)

add_executable( rtoz
  rtoz.cpp
)

target_link_libraries( rtoz
  rs
  cmdparser
  ITKCommon
  ITKIO
)

add_executable( randomizegraph
  randomizegraph.cpp
)

target_link_libraries( randomizegraph
  cmdparser
  ITKCommon
  ITKIO
  ITKStatistics
)

add_executable( fcmap
  fcmap.cpp
)

target_link_libraries( fcmap
  rs
  cmdparser
  ITKCommon
  ITKIO
)

add_executable( rmaptoz
  rmaptoz.cpp
)

target_link_libraries( rmaptoz
  rs
  cmdparser
  ITKCommon
  ITKIO
)

if( USE_BLAS )
add_executable( strengths strengths.cpp )
target_link_libraries( strengths flens blas cmdparser ITKCommon ITKIO )
endif( USE_BLAS )

if( USE_BLAS )
add_executable( clustercoefficient clustercoefficient.cpp )
target_link_libraries( clustercoefficient flens blas cmdparser ITKCommon ITKIO )
endif( USE_BLAS )

if( Boost_FOUND )
add_executable( fwsp fwsp.cpp )
target_link_libraries( fwsp cmdparser ITKCommon ITKIO )
endif( Boost_FOUND )

if( Boost_FOUND )
add_executable( hierarchy hierarchy.cpp )
target_link_libraries( hierarchy cmdparser flens blas ITKCommon ITKIO )
endif( Boost_FOUND )

