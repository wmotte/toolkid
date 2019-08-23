
#ifndef SIMPLEKMEANS_H_INCLUDED
#define SIMPLEKMEANS_H_INCLUDED

#include "DenseVector.h"

namespace gsis{

class SimpleKmeans{
    private:
        const int K; //number of clusters
        const int N; // number of points to be clustered
        const int dim; //dimension of space vector of points
        int e; // number of iterations
        double threshold;

        DenseVector ** points; // points to be clustered (a dynamic array of pointers to DenseVector)
        DenseVector ** muy; //cluster centers ( a dynamic array of pointers to DenseVector)
        int * assignment; // cluster assignment

    public:
        SimpleKmeans(int _K, int _N, int _dim, DenseVector ** _points); //constructor

        ~SimpleKmeans(); //destructor

        int * getAssignment();

        void clustering();

    private:

        void assignmentStep();

        double updateStep();

        //-----------------------------
        // IO method
        //-----------------------------
    public:
        /** input data: text file
         N dim
         Point 1: vector of dim dimension
         Point 2: vector of dim dimension
         .....
        Point N: vector of dim dimension*/

        static DenseVector ** readInputData(char * filename, int & N, int & dim);

        /** ouput cluster file: text file
        K dim: number of clusters and number of dimension
        <cluster 1>: vector of the cluster center
        <cluster 2>: vector of the cluster center
        ...
        <cluster K>: vector of the cluster center */

        void writeClusters(char * filename);

        /** ouput result file: text file
        N: number of points
        <cluster of point 1>
        ....
        <cluster of point N> */

        void writeResults(char * filename); // cluster id of the points


        /** output data:
         N dim
         Point 1: vector of cluster center to which this point belong
         Point 2: vector of cluster center to which this point belong
         .....
         Point N: vector of cluster center to which this point belong*/

        void writeClusterResult(char * filename); //replace points by their cluster muy
};

}

#endif // SIMPLEKMEANS_H_INCLUDED
