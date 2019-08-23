
#ifndef DENSEVECTOR_H_INCLUDED
#define DENSEVECTOR_H_INCLUDED

#include "Vector.h"

namespace gsis{

class DenseVector:public Vector{
    private:
        double * data;

    public:
        //-----------------------------
        // Constructor and Destructor
        //-----------------------------
        DenseVector(int size):Vector(size){
            data = new double[size];
        };

        virtual ~DenseVector(){
            delete [] data;
        }

        //-----------------------------
        // Other member methods
        //-----------------------------
        virtual double getElementAt(int index) const;

        virtual void setElementAt(int index, double val);

        virtual void addVector(const Vector & v);

        virtual Vector & subVector(const Vector & v);

        virtual Vector & divide(double val);

       virtual void copyFrom(const Vector & p);

        virtual void printTo(ostream & os);

        virtual void readFrom(istream & is);

        virtual double norm();

        virtual double distance(const Vector & v2);
};


}

#endif // DENSEVECTOR_H_INCLUDED
