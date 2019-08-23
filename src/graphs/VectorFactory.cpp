#include "VectorFactory.h"
#include "DenseVector.h"

#include <ctime>
#include <cstdlib>

using namespace gsis;

Vector * VectorFactory::zeros(int size){
    Vector * ret = new DenseVector(size);

    for (int i = 0; i < size; ++i){
        ret->setElementAt(i,0);
    }

    return ret;
}

Vector * VectorFactory::ones(int size){
    Vector * ret = new DenseVector(size);

    for (int i = 0; i < size; ++i){
        ret->setElementAt(i, 1);
    }

    return ret;
}

Vector * VectorFactory::uniform(int size){
    Vector * ret = new DenseVector(size);

    for (int i = 0; i < size; ++i){
        int temp = rand();

//        cout <<"temp:"<<temp << endl;

        double r = (double)temp/(double) RAND_MAX;

        ret->setElementAt(i, r);
    }
    return ret;
}
