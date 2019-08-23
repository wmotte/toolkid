#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <iomanip>
#include <new>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <string>

class Point {

public:

	Point( std::vector< double >& vs );

	~Point();

	double norm2DistanceSquared( Point p );

	std::vector< double >& getValues();

	unsigned int getDimension();

	std::string getStringValues();

	double getCoordinate( unsigned int c );

	void setCoordinate( unsigned int c, double value );

	Point operator+( Point point );

	Point operator+( double value );

	Point operator*( double value );

	Point operator/( double value );

private:
	std::vector< double > values;
	unsigned int dimension;
};

#endif

