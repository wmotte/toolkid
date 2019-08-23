
#include "SparseVector.h"
#include <map>
#include <math.h>
#include <string>
#include <cstdlib>
#include "Utils.h"

using namespace std;

using namespace gsis;

double SparseVector::norm(){
    double ret = 0;
    map<int, double>::iterator p = data->begin();

    for (; p != data->end(); ++p){
        ret += p->second * p->second;
    }

    return sqrt(ret);
}

double SparseVector::distance(const Vector & v2){
    //need to be improved
   double ret = 0.0;

    for (int i = 0; i < size; ++i){
        ret += pow(this->getElementAt(i) - v2.getElementAt(i),2);
    }

    return sqrt(ret);
}

double SparseVector::getElementAt(int index) const{
    map<int, double>::iterator p;

    p = data->find(index);
    if (p != data->end()){
        return p->second;
    }
    else return 0;
}

void SparseVector::addVector(const Vector & v2){
    //this needs to be improved
    for (int i = 0; i < size; ++i){
        double temp = getElementAt(i) + v2.getElementAt(i);
        setElementAt(i, temp);
    }
}


Vector & SparseVector::subVector(const Vector & v2){
    Vector * r = new SparseVector(size);

    for (int i = 0; i < size; ++i){
        double val = getElementAt(i)- v2.getElementAt(i);
        r->setElementAt(i,val);
    }

    return *r;
}

void SparseVector::setElementAt(int index, double val){
    map<int, double>::iterator p;
    p = data->find(index);

    if (val == 0){
        if (p != data->end())
            data->erase(p);
    }
    else {
        if (p != data->end())
            p->second = val; // update
        else
            data->insert(pair<int,double>(index,val)); //insert a new one
    }
}

void SparseVector::printTo(ostream & out){
    map<int, double>::iterator p = data->begin();

    for (; p != data->end(); ++p){
        out<< p->first<<":" << p->second << " ";
    }

}

void SparseVector::readFrom(istream & in){
    char ch;

    Utils::skipNewLineCharacter(in);
    while (true){
        Utils::moveToNextCharacter(in);

        ch = in.get();
        if (ch != '\n')
            in.putback(ch);
        else {
            break;
        }

        string temp;
        in >> temp;

        int loc = temp.find_first_of(":");
        int key = atoi(temp.substr(0, loc).c_str());
        double val = atof(temp.substr(loc + 1).c_str());

        setElementAt(key, val);
    }

}


Vector & SparseVector::divide(double val){
    //can be improved
    for (int i = 0; i < size; ++i){
        setElementAt(i, getElementAt(i)/val);
    }
    return *this;
}

void  SparseVector::copyFrom(const Vector & a){
    //can be improved
    if (a.getSize() != getSize()){
        cout << "Can not be assigned to a vector of different size";
    }
    else {
        for (int i = 0; i < size; ++i) //copying
            this->setElementAt(i, a.getElementAt(i));
    }
}



