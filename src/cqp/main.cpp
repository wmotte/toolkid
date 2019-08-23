#include "ConstrainedQuadraticProgramming.h"

#include <iostream>

/**
 * Test CQP.
 */
int main( int argc, char * const argv[] )
{
	cqp::MatrixType Q(2,2);
	Q(0,0) = 1.;
	Q(1,1) = 1.;

	cqp::VectorType q(2);
	q(0) = -2.;
	q(1) = -1;

	cqp::MatrixType A(1,2);
	A(0,0) = 3.;
	A(0,1) = 1.;


	cqp::VectorType b(1, 1.8 );

	cqp::MatrixType C(2,2);
	C(0,0) = 1.;
	C(1,1) = 1.;

	cqp::VectorType d(2, 0);

	unsigned int max_iter = 1000;

	cqp::ConstrainedQuadraticProgramming c = cqp::ConstrainedQuadraticProgramming( Q, q, A, b, C, d, max_iter );

	cqp::VectorType x = c.Solve();

	std::cout << "A: " << A << std::endl;
	std::cout << "x: " << x << std::endl;
	std::cout << "b: " << b << std::endl;

	return EXIT_SUCCESS;
}
