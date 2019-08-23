#include "Point.h"

Point::Point( int dim, double* p )
{
	dimension = dim;
	point = p;

}

Point* Point::AllocPoint( int dim, double* p )
{
	dimension = dim;
	point = new double[dim];

	for ( int i = 0; i < dim; i++ )
		point[i] = p[i];

	return this;
}

Point::Point( const Point& to_copy )
{
	dimension = to_copy.dimension;

	point = new double[dimension];
	if ( !point )
	{
		cerr << "Memory Allocation Failed in Point::Point" << endl;
		cerr << "Exitting the program..." << endl;
		exit( 1 );
	}
	for ( int i = 0; i < dimension; i++ )
		point[i] = to_copy.point[i];

}

Point::Point()
{
	dimension = 0;
	point = NULL;

}

Point::~Point()
{
	if ( point != NULL )
		delete[] point;

}

void Point::setPoint( double* p )
{
	point = p;
}

double Point::Norm2DistanceSquared( Point* p ) const
{
	double distance = 0;

	if ( dimension != p->getDimension() )
	{
		cerr << "Points are not of the same dimension, cannot calculate their norm-2 distance" << endl;
		return -1;
	}

	//calculate the Euclidean distance between this point and the passed point.
	for ( int i = 0; i < dimension; i++ )
	{
		distance += pow( ( point[i] - p->point[i] ), 2 );
	}

	return distance;
}

void Point::print( ostream* out ) const
{
	*out << "( ";
	for ( int i = 0; i < dimension; i++ )
	{
		if ( i != 0 )
			*out << ",";
		*out << setw( 7 ) << point[i];

	}
	*out << ")\n";

}

int Point::getDimension() const
{
	return dimension;
}

double* Point::getPoint() const
{
	return point;
}

int Point::operator==( const Point& p ) const
{

	for ( int i = 0; i < dimension; i++ )
	{
		if ( point[i] != p.point[i] )
			return 0;

	}

	return 1;
}

Point& Point::operator+( const Point& to_add ) const
{
	Point* result;

	if ( dimension != to_add.dimension )
	{
		cout << "can't add points of different dimension:" << endl;
		this->print( &cout );
		to_add.print( &cout );
		exit( 1 );
	}
	try
	{
		double* r = new double[dimension];

		for ( int i = 0; i < dimension; i++ )
		{
			r[i] = point[i] + to_add.point[i];
		}

		result = new Point( dimension, r );
	}
	catch ( int e )
	{
		cout << "Exception occured in Point::operator+: " << e << endl;
		cout << "Exiting program..." << endl;
		exit( 1 );
	}

	return *result;
}

Point& Point::operator*( double to_multiply ) const
{
	Point* result;

	try
	{
		double* r = new double[dimension];

		for ( int i = 0; i < dimension; i++ )
		{
			r[i] = point[i] * to_multiply;
		}

		result = new Point( dimension, r );

	} // try
	catch ( int e )
	{
		cout << "Exception occurred in Point::operator*: " << e << endl;
		cout << "Exiting program..." << endl;
		exit( 1 );
	}

	return *result;
}

Point& Point::operator/( double to_devide ) const
{

	if ( to_devide == 0 )
	{
		cout << " cannot devide a point by zero" << endl;
		exit( 0 );
	}

	double* r = new double[dimension];

	for ( int i = 0; i < dimension; i++ )
	{
		r[i] = point[i] / to_devide;
	}

	Point* result = new Point( dimension, r );

	return *result;
}

Point& Point::operator+( double to_add ) const
{
	double* r = new double[dimension];

	for ( int i = 0; i < dimension; i++ )
	{
		r[i] = point[i] + to_add;
	}

	Point* result = new Point( dimension, r );

	return *result;
}

Point& Point::operator-( double to_subtract ) const
{
	double* r = new double[dimension];

	for ( int i = 0; i < dimension; i++ )
	{
		r[i] = point[i] - to_subtract;
	}

	Point* result = new Point( dimension, r );

	return *result;
}

Point& Point::operator=( const Point& to_assign )
{
	int i;
	if ( this != &to_assign )
	{
		if ( dimension != to_assign.dimension )
		{
			dimension = to_assign.dimension;
			delete[] point;
			point = new double[dimension];
		}

		for ( i = 0; i < dimension; i++ )
			point[i] = to_assign.point[i];
	}

	return *this;
}

float Point::getCoordinate( int c )
{
	return point[c];
}

void Point::setCoordinate( int c, double value )
{
	point[c] = value;
}
