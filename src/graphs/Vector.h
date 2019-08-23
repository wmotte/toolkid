#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED

#include <iostream>

using namespace std;

namespace gsis{

class Vector{
    protected:
        const int size;

    public:
        //-----------------------------
        // Constructor and Destructor
        //-----------------------------

        Vector(int _size):size(_size){};

       //-----------------------------
        // Other member methods
        //-----------------------------
        virtual double getElementAt(int index) const = 0;

        virtual void setElementAt(int index, double val) = 0;

        virtual void addVector(const Vector & v) = 0;

        virtual Vector & subVector(const Vector & v) = 0;

        virtual Vector & divide(double val) = 0;

        virtual void printTo(ostream & os) = 0;

        virtual void readFrom(istream & is) = 0;

        int getSize()const {return size;}

        virtual void copyFrom(const Vector & p) = 0;

        virtual double norm() = 0;

        virtual double distance(const Vector & v2) = 0;

};

}
#endif // VECTOR_H_INCLUDED
