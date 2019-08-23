
#include <iostream>
#include <vector>
#include <boost/detail/algorithm.hpp>

/**
 * Test working boost header Xcode (17-08-2010).
 */
int main( int argc, char ** argv )
{

	std::vector< float > V( 5 );

	boost::iota( V.begin(), V.end(), 1 );

	/*
	 * Should print:
	 	 Value: 1
		 Value: 2
		 Value: 3
		 Value: 4
		 Value: 5
	 */
	for( unsigned int i = 0; i < V.size(); i++ )
		std::cout << "Value: " << V.at( i ) << std::endl;

	return 0;
}


