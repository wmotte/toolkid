

#include "SimpleKmeans.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "VectorFactory.h"
#include <cmath>

using namespace gsis;
using namespace std;

SimpleKmeans::SimpleKmeans(int _K, int _N, int _dim, DenseVector ** _points):K(_K),N(_N), dim(_dim){
    e = 100;
    threshold = 0.001;
    points = _points;

    //init cluster centers
    muy = new DenseVector*[K];
    for (int i = 0; i < K; ++i){
       //step 1: choose one center randomly among the data points
        muy[i] = new DenseVector(dim);
        double ran = (double) rand()/(double) RAND_MAX;
        int index = floor(ran * N);
        muy[i]->copyFrom(*points[index]);

    }

    assignment = new int[N];
}

SimpleKmeans::~SimpleKmeans(){
    for (int i = 0; i < K; ++i)
        delete muy[i];
    delete [] muy;

    delete [] assignment;
}

int * SimpleKmeans::getAssignment(){
    return assignment;
}

void SimpleKmeans::assignmentStep(){
    for (int i = 0; i < N; ++i){
        DenseVector * x = points[i];

        double minDistance = -1;
        for (int j = 0; j < K; ++j){
            DenseVector * c = muy[j];
            double distance = x->distance(*c);

            if (minDistance == -1){//first cluster
                minDistance = distance;
                assignment[i] = j;
            }
            else {
                if (distance < minDistance){
                    minDistance = distance;
                    assignment[i] = j;
                }
            }
        }
    }
}

double SimpleKmeans::updateStep(){
    //update cluster centers
    int * count = new int[K]; //count the number of elements in clusters

    DenseVector *m[K];
    for (int i = 0; i < K; ++i){
        m[i] = (DenseVector * ) VectorFactory::zeros(dim);
        count[i] = 0;
    }

    //counting
    for (int i = 0; i < N; ++i){
        count[assignment[i]] += 1;
        m[assignment[i]]->addVector(*points[i]);
       // muy[assignment[i]]->printTo(cout);
    }

    //calculate new cluster centers
    for (int i = 0; i < K; ++i){
        m[i]->divide(count[i]);
    }

    double maxShift = 0;
    for (int i = 0; i < K; ++i){
        double curShift = muy[i]->distance(*m[i]);
        if (curShift > maxShift)
            maxShift = curShift;

        muy[i]->copyFrom(*m[i]);
        delete m[i];
    }
    delete [] count;

    return maxShift;
}

void SimpleKmeans::clustering(){
    for (int iter = 0; iter < e; ++iter){
        cout << "Clustering step " << iter << ":";
        assignmentStep();
        double shift = updateStep();

        cout << shift << " shift" << endl;
        if (shift < threshold) break;
    }
}

//-----------------------------------
// IO methods
//-----------------------------------
DenseVector ** SimpleKmeans::readInputData(char * filename, int & N, int & dim){
    fstream fin(filename, fstream::in);
    char buffer[20];

    // read M, N and dim
    fin >> buffer;
    N = atoi(buffer);

    fin >> buffer;
    dim = atoi(buffer);

    DenseVector ** x = new DenseVector*[N];
    for (int i = 0; i < N; ++i){
        x[i] = new DenseVector(dim);
        for (int j = 0; j < dim; ++j){
            fin >> buffer;
            double temp = atof(buffer);
            x[i]->setElementAt(j,temp);
        }
    }

    fin.close();
    return x;
}

void SimpleKmeans::writeClusters(char * filename){
    fstream fout(filename, fstream::out);

    fout << N << endl;
    for (int i = 0; i < N; ++i){
        int idx = assignment[i];
        //c[i]->printTo(fout);
        fout << idx << endl;
    }

    fout.close();
}

void SimpleKmeans::writeClusterResult(char * filename){
    fstream fout(filename, fstream::out);

    //write N
    fout << N << "\t" << dim;
    for (int i = 0; i < N; ++i){
        int idx = assignment[i];
        DenseVector * c = muy[idx];
        c->printTo(fout);
        fout << endl;
    }

    fout.close();
}

void SimpleKmeans::writeResults(char * filename){
    fstream fout(filename, fstream::out);

    //write N
    fout << N << "\t" << dim << endl;
    for (int i = 0; i < N; ++i){
        int idx = assignment[i];
        //c[i]->printTo(fout);
        fout << idx << endl;
    }

    fout.close();
}
