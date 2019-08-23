/* file exmaple.C

 This file is an example on how uquadprog can be used
 by invoking solve_quadprog() function

 In order to compile this example, ublas boost library must be installed
 on your system.
 
 The test problem is the following:
 
 Given:
 G =  2.1 0.0 1.0   g0^T = [6.0 1.0 1.0]
 1.5 2.2 0.0
 1.2 1.3 3.1
 Solve:
 min f(x) = 1/2 x G x + g0 x
 s.t.
 x_1 + 2*x_2 + x_3 = -4

 x_1 >= 0
 x_2 >= 0
 x_3 >= 0
 -x_1 - x_2 >= -10
 
 The solution is x^T = [0 2 0] and f(x) = 6.4
 
 Author: Angelo Furfaro
 DEIS - University of Calabria, Italy
 a.furfaro@deis.unical.it
 http://www.lis.deis.unical.it/~furfaro
 
 LICENSE
 
 Copyright (2008) Angelo Furfaro
 


 This file is a part of uquadprog.

 uquadprog is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 uquadprog is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with uquadprog; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "uquadprog.hpp"
#include "vnlquadprog.hpp"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include <iostream>

#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>


/**
 * Ublas version.
 */
void ublasVersion()
{

	double rd = boost::math::ellint_rd( 1.2, 1.3, 1.4 );
	double rf = boost::math::ellint_rf( 1.2, 1.3, 1.4 );

	std::cout << "Rd: " << rd << std::endl;
	std::cout << "Rf: " << rf << std::endl;

	std::cout << "ublas version" << std::endl;
	matrix< double > G( 3, 3 );
	vector< double > g0( 3 );
	vector< double > x( 3 );

	G( 0, 0 ) = 2.1;
	G( 0, 1 ) = 0.0;
	G( 0, 2 ) = 1.2; // 1.0;

	G( 1, 0 ) = 0.0; // 1.5;
	G( 1, 1 ) = 2.2;
	G( 1, 2 ) = 1.3; // 0.0;

	G( 2, 0 ) = 1.2;
	G( 2, 1 ) = 1.3;
	G( 2, 2 ) = 3.1;

	g0( 0 ) = 6.0;
	g0( 1 ) = 1.0;
	g0( 2 ) = 1.0;

	// equality constraints

	matrix< double > CE( 3, 1 );
	vector< double > ce0( 1 );

	CE( 0, 0 ) = 1.0;
	CE( 1, 0 ) = 2.0;
	CE( 2, 0 ) = -1.0;

	ce0( 0 ) = -4;

	// inequality constraints.
	matrix< double > CI( 3, 4 );
	vector< double > ci0( 4 );

	CI( 0, 0 ) = 1.0;
	CI( 0, 1 ) = 0.0;
	CI( 0, 2 ) = 0.0;
	CI( 0, 3 ) = -1.0;
	CI( 1, 0 ) = 0.0;
	CI( 1, 1 ) = 1.0;
	CI( 1, 2 ) = 0.0;
	CI( 1, 3 ) = -1.0;
	CI( 2, 0 ) = 0.0;
	CI( 2, 1 ) = 0.0;
	CI( 2, 2 ) = 1.0;
	CI( 2, 3 ) = 0.0;

	ci0( 0 ) = 0.0;
	ci0( 1 ) = 0.0;
	ci0( 2 ) = 0.0;
	ci0( 3 ) = 10.0;

	std::cout << "f: " << solve_quadprog( G, g0, CE, ce0, CI, ci0, x ) << std::endl;
	std::cout << "x: ";
	for ( unsigned int i = 0; i < x.size(); i++ )
		std::cout << x( i ) << ' ';
	std::cout << std::endl;
}

/**
 * vnl version.
 */
void vnlVersion()
{
	std::cout << "vnl version" << std::endl;

	typedef vnl_matrix< double > MatrixType;
	typedef vnl_vector< double > VectorType;

	MatrixType G( 3, 3 );
	VectorType g0( 3 );
	VectorType x( 3 );

	G( 0, 0 ) = 2.1;
	G( 0, 1 ) = 0.0;
	G( 0, 2 ) = 1.2; // 1.0;

	G( 1, 0 ) = 0.0; // 1.5;
	G( 1, 1 ) = 2.2;
	G( 1, 2 ) = 1.3; // 0.0;

	G( 2, 0 ) = 1.2;
	G( 2, 1 ) = 1.3;
	G( 2, 2 ) = 3.1;

	g0( 0 ) = 6.0;
	g0( 1 ) = 1.0;
	g0( 2 ) = 1.0;

	// equality constraints

	MatrixType CE( 3, 1 );
	VectorType ce0( 1 );

	CE( 0, 0 ) = 1.0;
	CE( 1, 0 ) = 2.0;
	CE( 2, 0 ) = -1.0;

	ce0( 0 ) = -4;

	// inequality constraints.
	MatrixType CI( 3, 4 );
	VectorType ci0( 4 );

	CI( 0, 0 ) = 1.0;
	CI( 0, 1 ) = 0.0;
	CI( 0, 2 ) = 0.0;
	CI( 0, 3 ) = -1.0;
	CI( 1, 0 ) = 0.0;
	CI( 1, 1 ) = 1.0;
	CI( 1, 2 ) = 0.0;
	CI( 1, 3 ) = -1.0;
	CI( 2, 0 ) = 0.0;
	CI( 2, 1 ) = 0.0;
	CI( 2, 2 ) = 1.0;
	CI( 2, 3 ) = 0.0;

	ci0( 0 ) = 0.0;
	ci0( 1 ) = 0.0;
	ci0( 2 ) = 0.0;
	ci0( 3 ) = 10.0;


	std::cout << "f: " <<  qp::solve_quadprog( G, g0, CE, ce0, CI, ci0, x ) << std::endl;
	std::cout << "x: ";
	for ( unsigned int i = 0; i < x.size(); i++ )
		std::cout << x( i ) << ' ';
	std::cout << std::endl;

}

/**
 * main.
 */
int main( int argc, char** argv )
{
	ublasVersion();
	vnlVersion();
}



