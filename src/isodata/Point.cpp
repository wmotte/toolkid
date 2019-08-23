#include "Point.h"
#include <sstream>

Point::Point( std::vector< double >& vs ) {
	values = vs;
	dimension = values.size();
}

Point::~Point() {
}

/**
 * Return Euclidean distance squared.
 */
double Point::norm2DistanceSquared( Point p ) {

	std::vector< double > new_values = p.getValues();

	double distance = 0;

	if ( dimension != p.getDimension() ) {
		std::cerr << "Points are not of the same dimension!" << std::endl;
		return -1;
	}

	// Euclidean distance...
	for ( unsigned int i = 0; i < dimension; i++ ) {
		distance += pow( ( values[i] - new_values[i] ), 2 );
	}
	return distance;
}

/**
 * Get values.
 */
std::vector< double >& Point::getValues() {
	return values;
}

/**
 * Return dimension.
 */
unsigned int Point::getDimension() {
	return dimension;
}

/**
 * Print all values for point.
 */
std::string Point::getStringValues() {
	std::stringstream ss;
	ss << "[ ";

	for ( unsigned int i = 0; i < dimension; i++ ) {
		ss << values[i] << " ";
	}

	ss << "]" << std::endl;
	return ss.str();
}

double Point::getCoordinate( unsigned int c ) {
	return values[c];
}

void Point::setCoordinate( unsigned int c, double value ) {
	values[c] = value;
}

Point Point::operator+( Point point ) {
	std::vector< double > new_values = point.getValues();
	for ( unsigned int i = 0; i < values.size(); i++ ) {
		values[i] += new_values[i];
	}
	return Point( values );
}

Point Point::operator+( double value ) {
	for ( unsigned int i = 0; i < values.size(); i++ ) {
		values[i] += value;
	}
}

Point Point::operator*( double value ) {
	for ( unsigned int i = 0; i < values.size(); i++ ) {
		values[i] *= value;
	}
}

Point Point::operator/( double value ) {
	for ( unsigned int i = 0; i < values.size(); i++ ) {
		values[i] /= value;
	}
}

