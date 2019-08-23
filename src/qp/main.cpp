/*
 File main.cc
 
 This file contains just an example on how to set-up the matrices for using with
 the solve_quadprog() function.
 
 The test problem is the following:
 
 Given:
 H =  4 -2   f' = [6 0]
 -2  4
 
 Solve:
 min f(x) = 1/2 x H x + f x

 s.t.
 (b) => 	x_1 + x_2 = 3

 (d) => 	x_1 >= 0
 	 	 	x_2 >= 0
 	 	 	x_1 + x_2 >= 2
 
 The solution is x' = [1 2] and f(x) = 12

 */

#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"

using namespace QuadProgPP;

/**
 * Test QP.
 */
int main( int argc, char * const argv[] )
{
	Matrix< double > H, A, C;
	Vector< double > f, b, d, x;
	int n, m, p;

	char ch;

	n = 2;
	H.resize( n, n );
	{
		std::istringstream is( "4, -2,"
			"-2, 4 " );

		for ( int i = 0; i < n; i++ )
			for ( int j = 0; j < n; j++ )
				is >> H[i][j] >> ch;
	}

	f.resize( n );
	{
		std::istringstream is( "6.0, 0.0 " );

		for ( int i = 0; i < n; i++ )
			is >> f[i] >> ch;
	}

	m = 1;
	A.resize( n, m );
	{
		std::istringstream is( "1.0, "
			"1.0 " );

		for ( int i = 0; i < n; i++ )
			for ( int j = 0; j < m; j++ )
				is >> A[i][j] >> ch;
	}

	b.resize( m );
	{
		std::istringstream is( "-3.0 " );

		for ( int j = 0; j < m; j++ )
			is >> b[j] >> ch;
	}

	p = 3;
	C.resize( n, p );
	{
		std::istringstream is( "1.0, 0.0, 1.0, "
			"0.0, 1.0, 1.0 " );

		for ( int i = 0; i < n; i++ )
			for ( int j = 0; j < p; j++ )
				is >> C[i][j] >> ch;
	}

	d.resize( p );
	{
		std::istringstream is( "0.0, 0.0, -2.0 " );

		for ( int j = 0; j < p; j++ )
			is >> d[j] >> ch;
	}
	x.resize( n );

	// ------------------------------------
	// Ax - b  = 0
	// Cx - d >= 0
	// ------------------------------------

	std::cout << std::endl << "H: " << std::endl;
	std::cout << H << std::endl;

	solve_quadprog( H, f, A, - b, C, d, x );

	std::cout << "b: " << std::endl;
	for( unsigned int i = 0; i < b.size(); i++ )
		std::cout << b[i] << std::endl;

	std::cout << std::endl << "d: " << std::endl;
	for( unsigned int i = 0; i < d.size(); i++ )
		std::cout << d[i] << std::endl;

	std::cout << std::endl << "f: " << std::endl;
		for( unsigned int i = 0; i < f.size(); i++ )
			std::cout << f[i] << std::endl;

	std::cout << std::endl << "A: " << std::endl;
	std::cout << A << std::endl;

	std::cout << std::endl << "C: " << std::endl;
	std::cout << C << std::endl;

	std::cout << std::endl << "x (solution): " << std::endl;
	for( unsigned int i = 0; i < x.size(); i++ )
		std::cout << x[i] << std::endl;





}
