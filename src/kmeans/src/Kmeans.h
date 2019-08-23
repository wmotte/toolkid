
// From: Arthur, D. and Vassilvitskii, S., "k-means++: the advantages of careful seeding",
// Proc. of ACM-SIAM, 2007


#ifndef KMEANS_H_INCLUDED
#define KMEANS_H_INCLUDED

#include "DenseVector.h"
#include "SparseVector.h"
#include <vector>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

namespace gsis{

	typedef vnl_vector< double > VectorType;
	typedef vnl_matrix< double > MatrixType;

class Kmeans{
    private:
        const int K;      // number of clusters
        const long N;     // number of points to be clustered
        const int dim;    //dimension of points and clusters
        int e;            //number of iterations

        double threshold;   //the stop condition for Kmeans
        bool kmeanspp;      //whether or not initialize Kmeans according Kmeans++

        DenseVector ** muy;     // cluster centers (size of K) - array of pointers
        double ** l;            //lower bound (size of N*K)
        Vector ** x;      // points to be clustered (size of N) - array of pointers
        double * u;             //upper bound of distance from x to cx(x) (size of N)
        double ** d;            // distances among cluster centers

        bool * r; //size of N
        long * assignment; //size of N, index of cluster to which a point is currently assigned

        double * s; //s(c) = min_c' d(c,c') - size of K

        void seedKMeansClusters();

        double initialize(int _K);

    public:
        //-------------------------------
        // Constructor & Destructor
        //-------------------------------
        Kmeans(int _K, long _N, int _dim, Vector ** _x); //constructor

        ~Kmeans(); //destructor

        //--------------------------------
        // method
        //--------------------------------
        void clustering();

        //get, set methods
        long * getAssignment();

        const int getK() const;

        const int getN()const;

        const int getDim() const;

        void setNIter(int _niter);

        void setThreshold(double _threshold);

        void setKmeansPP(bool _kmeanspp);

        //------------------------------------
        // IO Methods
        //------------------------------------
            /** input data: text file
         N dim
         Point 1: vector of dim dimension
         Point 2: vector of dim dimension
         .....
        Point N: vector of dim dimension*/

        static Vector ** readInputData(char * filename, int & N, int & dim, bool sparse = true);

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

        /**
         * Get vector with assigned clusters.
         */
        VectorType getResults();

        /**
         * Convert vnl MatrixType into Vector**.
         */
        static Vector ** readInputData( const MatrixType& matrix, int & N, int & dim );
};


}


#endif // KMEANS_H_INCLUDED
