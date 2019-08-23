#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <iomanip>
#include <new>
#include <cmath>
#include <cstdlib>

using namespace std;

class Point
{
public:

	Point();
	Point( int dim, double* p );
	Point* AllocPoint( int dim, double* p );
	Point( const Point& );
	~Point();
	double Norm2DistanceSquared( Point* p ) const;
	void print( ostream* out ) const;
	int getDimension() const;

	double* getPoint() const;
	void setPoint( double* p );

	Point& operator+( const Point& ) const;
	Point& operator+( double ) const;
	Point& operator-( double ) const;
	Point& operator*( double ) const;
	Point& operator/( double ) const;
	Point& operator=( const Point& );
	int operator==( const Point& ) const;
	float getCoordinate( int c );
	void setCoordinate( int c, double value );

private:

	double* point;
	int dimension;
};
#endif

