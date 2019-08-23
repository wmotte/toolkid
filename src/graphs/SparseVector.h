

#ifndef SPARSEVECTOR_H_INCLUDED
#define SPARSEVECTOR_H_INCLUDED
#include <map>
#include "Vector.h"

using namespace std;
namespace gsis{

class SparseVector:public Vector{
    private:
        map<int, double>  * data;

    public:
        //-----------------------------
        // Constructor and Destructor
        //-----------------------------
        SparseVector(int size):Vector(size){
            data = new map<int,double>();
        };

        virtual ~SparseVector(){
            delete data;
        }

         //-----------------------------
        // Other member methods
        //-----------------------------
        virtual double getElementAt(int index) const;

        virtual void setElementAt(int index, double val);

        virtual void addVector(const Vector & v2);

        virtual Vector & subVector(const Vector & v2);

        virtual Vector & divide(double val);

        virtual void copyFrom(const Vector & p);

        virtual void printTo(ostream & os);

        virtual void readFrom(istream & is);

        virtual double norm();

        virtual double distance(const Vector & v2);
};
}
#endif // SPARSEVECTOR_H_INCLUDED
