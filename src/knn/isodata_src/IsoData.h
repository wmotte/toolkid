#ifndef ISODATA_H
#define ISODATA_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <cfloat>
#include <limits>

#include "Point.h"
#include "KMrand.h"
#include "KM_ANN.h"
#include "KMfilterCenters.h"

// consts
static const float K = 0.5;
static const int MINCLUS = 1;

extern vector< vector< int > > clusters;
extern vector< KMpoint > kcCenters;

/**
 * PairDistanceNode contains the pair wise distance between two cluster centers c1, and c2.
 * We will need to save the distance as well as the clusters number of which c1, and c2 are
 * their centers.
 */
struct PairDistanceNode
{
	double dist;
	int c1;
	int c2;

	bool operator<( PairDistanceNode to_compare ) const
	{
		return ( dist < to_compare.dist );
	}

	bool operator<( double to_compare ) const
	{
		return ( dist < to_compare );
	}
};

/**
 * Modified from: "A Fast Implementation of the ISODATA Clustering Algorithm",
 * by 'Nargess Memarsadeghi, David M. Mount, Nathan S. Netanyahu, and Jacqueline Le Moigne'.
 *
 * [http://www.cs.umd.edu/~mount/Projects/ISODATA/]
 *
 * IsoData. @author: wim@invivonmr.uu.nl; Image Sciences Institute, UMC Utrecht.
 */
class IsoData
{

public:

	IsoData( int pixels, int nImages, int nClusters, int minNumberOfSamplesInCluster );
	~IsoData();
	void DeleteCenters();
	void BuildKMfilterCenters();
	void SetNumClusters();
	void SetFilterCenters();
	void SetImageCenters();
	void ReadImage( string* );
	void WriteClassifiedImage( string );
	void UpdateCenters();
	void SetPoints( const KMpointArray& all );
	int GetNumCenters();
	int SampleSize();
	bool WasDeleted();
	void SampleCenters();
	void SamplePoints( double );
	void PutInCluster();
	void PostAnalyzeClusters( bool verbose );
	void CalculateAverageDistances();
	double OverallAverageDistances();
	void CalculateSTDVector();
	void CalculateVmax();
	vector< int > ShouldSplit( double stdv );
	void Split( vector< int > to_split, bool verbose );
	void ComputeCenterDistances();
	vector< PairDistanceNode > FindLumpCandidates( double lump, int MAXPAIR );
	void Lump( const vector< PairDistanceNode >& to_lump, bool verbose );
	void PreFinalClustering();
	void GenerateReport( ostream* );
	std::vector< std::vector< int > > GetClusters();

private:

	int NumBands;
	int ImageSizeInByte;
	int NumClusters;
	unsigned int NumSamples;
	bool Deleted;
	Point** allPoints;
	vector< Point* > centers;
	KMdata* data;
	KMfilterCenters* filter;
	vector< int > centers_to_keep;
	vector< int > samples;
	double** STDVector;
	double* Vmax;
	int* Vmax_index;
	double* average_distances;
	double OverallD;
	vector< PairDistanceNode > CenterDistances;
	int FindCenter( Point* to_find );

};
#endif

