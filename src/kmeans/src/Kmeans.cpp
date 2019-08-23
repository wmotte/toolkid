#include "Kmeans.h"

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "VectorFactory.h"
#include <math.h>

using namespace gsis;
using namespace std;

Kmeans::Kmeans( int _K, long _N, int _dim, Vector **_x ) :
	K( _K ), N( _N ), dim( _dim )
{

	//initialize data
	cout << "Initialization " << endl;
	x = _x;

	d = new double *[K];
	for ( int i = 0; i < K; ++i )
		d[i] = new double[K];

	//  cout << "Initialization 2" << endl;
	//initialize l (dynamic 2-dim array)
	l = new double*[N];
	for ( int i = 0; i < N; ++i )
		l[i] = new double[K];

	// cout << "Initialization 3" << endl;
	//initialize u
	u = new double[N];

	//initialize r
	// cout << "Initialization 4" << endl;
	r = new bool[N];
	assignment = new long[N];
	e = 100;
	threshold = 2;
	kmeanspp = false;

	s = new double[K];

	//cout << "Initialization 5" << endl;
	//initialize cluster centers
	seedKMeansClusters();
}

//Seed Kmeans according to Kmeans++
void Kmeans::seedKMeansClusters()
{
	if ( kmeanspp )
	{//Seed Kmeans according to Kmeans++
		muy = new DenseVector*[K];
		for ( int i = 0; i < K; ++i )
			muy[i] = new DenseVector( dim );

		//step 1: choose one center randomly among the data points
		double ran = (double) rand() / (double) RAND_MAX;
		int index = floor( ran * N );
		muy[0]->copyFrom( *x[index] );

		// repeatedly choose more center
		for ( int _K = 1; _K < K; ++_K )
		{
			cout << "seeding " << _K << endl;
			//step 2: for each data point,
			//compute u(x): distance between x and the nearest center that has already been chosen
			double total_cost = initialize( _K );

			//step 3: add one new data point as new center
			double cut_off = (double) rand() / (double) RAND_MAX * total_cost;
			double cur_cost = 0;
			for ( int i = 0; i < N; ++i )
			{
				cur_cost += u[i];
				if ( cur_cost > cut_off )
				{
					muy[_K]->copyFrom( *x[i] );
					break;
				}
			}
		}
	} else
	{//otherwise
		muy = new DenseVector*[K];
		for ( int i = 0; i < K; ++i )
		{
			//step 1: choose one center randomly among the data points
			muy[i] = new DenseVector( dim );
			double ran = (double) rand() / (double) RAND_MAX;
			int index = floor( ran * N );
			muy[i]->copyFrom( *x[index] );
		}
	}

	initialize( K );
}

Kmeans::~Kmeans()
{
	delete[] assignment;
	delete[] r;
	delete[] u;
	delete[] s;

	for ( int i = 0; i < N; ++i )
	{
		delete[] l[i];
	}
	delete[] l;

	for ( int i = 0; i < K; ++i )
		delete muy[i];
	delete[] muy;

	for ( int i = 0; i < K; ++i )
		delete d[i];
	delete d;
}

double Kmeans::initialize( int _K )
{
	//set the lowerbound l[x][c] = 0 for each point x and center c
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < _K; ++j )
			l[i][j] = 0;
	}

	//calculate d(c,c') between any pairs of cluster centers c and c'
	for ( int i = 0; i < _K; ++i )
		for ( int j = i; j < _K; ++j )
		{
			if ( i == j )
			{
				d[i][j] = 0;
				continue;
			}

			double _d = muy[i]->distance( *muy[j] );
			d[i][j] = _d;
			d[j][i] = _d;
		}

	//assign each x to its closest initial center assignment= argmin_c d(x,c)
	// using lemma 1 to avoid redundant distance calculations
	double total_cost = 0;
	for ( int i = 0; i < N; ++i )
	{
		r[i] = true;
		assignment[i] = -1; // store the cluster with minimum distance
		double minD = -1;
		bool dis_cal = false; //whether we have to calculate the distance or not

		for ( int j = 0; j < _K; ++j )
		{
			if ( assignment[i] == -1 )
			{
				dis_cal = true;
			} else
			{
				//get the current assigned center
				int _c = assignment[i];
				double temp = 0.5 * d[j][_c];
				if ( temp >= minD )
					dis_cal = false;
				else
					dis_cal = true;
			}

			if ( dis_cal == true )
			{
				double _d = x[i]->distance( *muy[j] );
				l[i][j] = _d;

				if ( minD >= _d || minD == -1 )
				{
					assignment[i] = j;
					minD = _d;
				}
			}
		}

		u[i] = minD;
		total_cost += u[i];
	}

	return total_cost;
}

void Kmeans::clustering()
{
	//repeat until convergence
	for ( int iter = 0; iter < e; ++iter )
	{
		//cout << "clustering step: " << iter << endl;

		// for all centers c and c', compute d(c,c'), and s(c) = 1/2 min (c'#c)d(c,c')
		for ( int _c = 0; _c < K; ++_c )
		{
			s[_c] = -1;

			for ( int __c = 0; __c < K; ++__c )
			{
				if ( _c == __c )
					d[_c][__c] = 0;
				else
				{
					double temp = muy[_c]->distance( *muy[__c] );
					d[_c][__c] = temp;

					if ( s[_c] == -1 )
					{
						s[_c] = 0.5 * d[_c][__c];
					} else if ( 0.5 * d[_c][__c] < s[_c] )
					{
						s[_c] = 0.5 * d[_c][__c];
					}
				}
			}

		}

		// Identify all points such that u(x) <= s(c(x))
		// For all remaining points and centers c such that
		// c# c(x)
		// u(x) > l(x,c)
		// u(x) > 1/2 d(c(x),c)
		for ( int i = 0; i < N; ++i )
		{
			if ( u[i] <= s[assignment[i]] )
			{
				// do nothing, keep clusters
				//  cout << "do nothing with" << i << endl;
			} else
			{
				double dx_cx = 0;
				for ( int _c = 0; _c < K; _c++ )
				{
					// cout << l[i][_c] << endl;
					if ( _c != assignment[i] && u[i] > l[i][_c] && u[i] > 0.5 * d[assignment[i]][_c] )
					{
						if ( r[i] )
						{
							dx_cx = x[i]->distance( *muy[assignment[i]] );
							u[i] = dx_cx;
							r[i] = false;
						} else
							dx_cx = u[i];

						if ( dx_cx > l[i][_c] || dx_cx > 0.5 * d[assignment[i]][_c] )
						{
							//compute d(x,c)
							double temp = x[i]->distance( *muy[_c] );
							l[i][_c] = temp;
							if ( temp < dx_cx )
							{
								assignment[i] = _c;
								u[i] = temp;
							}
						}
					}
				}
			}
		}

		// for each center c, let m(c) be the mean of the points assigned to c
		DenseVector *m[K]; //initialize outside the main loop to speed up
		int count[K]; //initialize outside the main loop to speed up

		for ( int i = 0; i < K; ++i )
		{
			m[i] = (DenseVector *) VectorFactory::zeros( dim );
			count[i] = 0;
		}

		for ( int i = 0; i < N; ++i )
		{
			m[assignment[i]]->addVector( *x[i] );
			count[assignment[i]] += 1;
		}

		for ( int i = 0; i < K; ++i )
		{
			//cout << i << ":" << count[i] << endl;
			m[i]->divide( count[i] );
		}

		// for each point x and center c, update l(x,c)
		// calculate distances between muy and new center m of the same cluster in advance
		double centerShift[K]; // initialize outside the main loop to speed up
		double maxShift = 0;
		for ( int i = 0; i < K; ++i )
		{
			centerShift[i] = muy[i]->distance( *m[i] );
			if ( centerShift[i] > maxShift )
				maxShift = centerShift[i];
		}

		for ( int i = 0; i < N; ++i )
		{
			for ( int j = 0; j < K; ++j )
			{
				double temp = l[i][j] - centerShift[j];
				if ( temp > 0 )
					l[i][j] = temp;
				else
					l[i][j] = 0;
			}

			// for each point
			u[i] = u[i] + centerShift[assignment[i]];
			r[i] = true;
		}

		//replace each center c by m(c)
		for ( int i = 0; i < K; ++i )
		{
			muy[i]->copyFrom( *m[i] );

			delete m[i];
		}
		//cout << "max shift:" << maxShift << endl;
		if ( maxShift < threshold )
			break;

	}//end of the main loop
}

long * Kmeans::getAssignment()
{
	return assignment;
}

const int Kmeans::getK() const
{
	return K;
}

const int Kmeans::getN() const
{
	return N;
}

const int Kmeans::getDim() const
{
	return dim;
}

void Kmeans::setNIter( int _niter )
{
	e = _niter;
}

void Kmeans::setKmeansPP( bool _kmeanspp )
{
	kmeanspp = _kmeanspp;
}

void Kmeans::setThreshold( double _threshold )
{
	threshold = _threshold;
}

//------------------------------------
// IO methods
//------------------------------------
Vector ** Kmeans::readInputData( char * filename, int & N, int & dim, bool sparse )
{
	//cout << "reading data from " << filename << endl;
	fstream fin( filename, fstream::in );
	char buffer[100];

	// read M, N and dim
	fin >> buffer;
	N = atoi( buffer );
	fin >> buffer;
	dim = atoi( buffer );

	//cout << "dim " << N << " " << dim << endl;
	Vector ** x = new Vector*[N];
	if ( sparse )
	{
		for ( int i = 0; i < N; ++i )
		{
			//cout << "point " << i << endl;
			x[i] = new SparseVector( dim );
			/*for (int j = 0; j < dim; ++j){
			 fin >> buffer;
			 double temp = atof(buffer);
			 x[i]->setElementAt(j,temp);
			 }*/
			x[i]->readFrom( fin );
			// x[i]->printTo(cout);
		}
		//cout << "finish reading data";
	}
	else
	{
		for ( int i = 0; i < N; ++i )
		{
			x[i] = new DenseVector( dim );
			/*for (int j = 0; j < dim; ++j){
			 fin >> buffer;
			 double temp = atof(buffer);
			 x[i]->setElementAt(j,temp);
			 }*/
			x[i]->readFrom( fin );
		}
	}

	fin.close();
	return x;
}

/**
 * Read double matrix into Vector **.
 */
Vector ** Kmeans::readInputData( const MatrixType& matrix, int& N, int& dim )
{
	if( matrix.empty() )
	{
		std::cerr << "*** ERROR ***: matrix is empty (Kmeans::readInputData)!" << std::endl;
		exit( EXIT_FAILURE );
	}
	//N = matrix.rows();
	//dim = matrix.cols();
	Vector ** x = new Vector*[N];

	for( int i = 0; i < N; ++i )
	{
		x[i] = new DenseVector( dim );
		for( int j = 0; j < dim; ++j )
		{
			x[i]->setElementAt( j, matrix[i][j] );
		}
	}
	return x;
}

void Kmeans::writeResults( char * filename )
{
	fstream fout( filename, fstream::out );

	//write N
	fout << N << "\t" << dim << endl;
	for ( int i = 0; i < N; ++i )
	{
		int idx = assignment[i];
		fout << idx << endl;
	}

	fout.close();
}

/**
 * Return vector with assigned clusters to points.
 */
VectorType Kmeans::getResults()
{
	VectorType result( N );

	for ( int i = 0; i < N; ++i )
		result[i] = assignment[i] + 1; // no zero-based clusters...

	return result;
}

void Kmeans::writeClusters( char * filename )
{
	fstream fout( filename, fstream::out );

	//write K
	fout << K << endl;
	for ( int i = 0; i < K; ++i )
	{
		muy[i]->printTo( fout );
		fout << endl;
	}
	fout.close();
}

void Kmeans::writeClusterResult( char * filename )
{
	fstream fout( filename, fstream::out );

	//write N
	fout << N << "\t" << dim << endl;
	for ( int i = 0; i < N; ++i )
	{
		int idx = assignment[i];
		DenseVector * c = muy[idx];
		c->printTo( fout );
		fout << endl;
	}

	fout.close();
}

