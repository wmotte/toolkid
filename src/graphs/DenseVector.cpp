
#include "DenseVector.h"
#include <math.h>
#include <iostream>

using namespace std;
using namespace gsis;

double DenseVector::getElementAt(int index) const{
    return data[index];
}

void DenseVector::setElementAt(int index, double val){
    data[index] = val;
}

void DenseVector::addVector(const Vector & v2){
     for (int i = 0; i < size; ++i){
        //cout << "v2:" << v2.getElementAt(i);
        double temp = data[i] + v2.getElementAt(i);
        setElementAt(i, temp);
    }
}


Vector & DenseVector::subVector(const Vector & v2){
    DenseVector * r = new DenseVector(size);

    for (int i = 0; i < size; ++i){
        double val = data[i] - v2.getElementAt(i);
        r->setElementAt(i,val);
    }

    return *r;
}

Vector & DenseVector::divide(double val){
    for (int i = 0; i < size; ++i){
        setElementAt(i, getElementAt(i)/val);
    }
    return *this;
}

void DenseVector::copyFrom(const Vector & a){
    if (a.getSize() != getSize()){
        cout << "Can not be assigned to a vector of different size:" << a.getSize() << ":" << getSize();
    }
    else {
        for (int i = 0; i < size; ++i) //copying
            this->setElementAt(i, a.getElementAt(i));
    }
}

double DenseVector::norm(){
    double ret = 0.0;

    for (int i = 0; i < size; ++i){
        ret += data[i] * data[i];
    }
    return sqrt(ret);
}

double DenseVector::distance(const Vector & v2){
    double ret = 0.0;

    for (int i = 0; i < size; ++i){
        ret += pow(data[i] - v2.getElementAt(i),2);
    }

    return sqrt(ret);
}

void DenseVector::printTo(ostream & out){
    for (int i = 0; i < size; ++i){
        out << data[i]<<"\t";
    }

    //out << "\n";
}

void DenseVector::readFrom(istream & in){
    double val;

    for (int i = 0; i < size; ++i){
        in >> val;

        data[i] = val;
    }
}
