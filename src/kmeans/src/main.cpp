
#include <iostream>
#include "Kmeans.h"
#include "SimpleKmeans.h"
#include "SparseVector.h"
#include "DenseVector.h"
#include "VectorFactory.h"
#include <string>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <math.h>

//using namespace gsis;

/**
 * Test.
 */
int main( int argc, char * argv[] )
{
	if ( argc != 4 )
	{
		std::cout << "kmeans++2" << std::endl;
		exit( 1 );
	}

	srand( (unsigned) time( 0 ) );

	gsis::MatrixType matrix( 4,2 );

	matrix[0][0] = 1.2;
	matrix[0][1] = 1.3;
	matrix[0][2] = 4.2;
	matrix[0][3] = 5.2;

	matrix[1][0] = 7.2;
	matrix[1][1] = 8.3;
	matrix[1][2] = 6.2;
	matrix[1][3] = 0.2;

	int N = matrix.rows();
	int dim = matrix.cols();
	gsis::Vector ** x = gsis::Kmeans::readInputData( matrix, N, dim );
	int K = 2;
	double threshold = exp( -10 );

	gsis::Kmeans cluster( K, N, dim, x );
	cluster.setThreshold( threshold );
	cluster.setKmeansPP( true ); // use kmeans++ seeding
	cluster.clustering();

	gsis::VectorType assigned = cluster.getResults();
	std::cout << "N: " << assigned.size() << std::endl;
	std::cout << "Threshold: " << threshold << std::endl;
	std::cout << "K: " << K << std::endl;
	for( unsigned int i = 0; i < assigned.size(); i++ )
	{
		std::cout << "Cluster: " << assigned[i] << std::endl;
	}

	return 0;
}
