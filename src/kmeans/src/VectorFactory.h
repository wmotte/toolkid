#ifndef VECTORFACTORY_H_INCLUDED
#define VECTORFACTORY_H_INCLUDED

#include "Vector.h"

namespace gsis{

class VectorFactory{
public:
    static Vector * zeros(int size);

    static Vector * ones (int size);

    static Vector * uniform(int size);
};

}

#endif // VECTORFACTORY_H_INCLUDED
