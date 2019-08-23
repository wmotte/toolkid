#include "IsoData.h"

int sample_seed = 0;

/**
 * Grab clusters.
 */
std::vector< std::vector< int > > IsoData::GetClusters()
{
	return clusters;
}

/**
 * Constructor IsoData.
 */
IsoData::IsoData( int pixels, int nImages, int nClusters, int minNumberOfSamplesInCluster )
{
	sample_seed = time( 0 );

	NumBands = nImages;
	NumClusters = nClusters;
	NumSamples = minNumberOfSamplesInCluster;
	Deleted = false;
	ImageSizeInByte = pixels;
	STDVector = NULL;
	Vmax = NULL;
	Vmax_index = NULL;

	allPoints = new Point*[ImageSizeInByte];

	if ( !allPoints )
	{
		cerr << "Error 3: Memory Allocation Failed in IsoData::Image" << endl;
		exit( 1 );
	}

	for ( int i = 0; i < ImageSizeInByte; i++ )
	{
		allPoints[i] = NULL;
	}
}

/**
 * Destructor.
 */
IsoData::~IsoData()
{
	DeleteCenters();

	delete[] allPoints;
	delete[] average_distances;

	allPoints = NULL;
	samples.clear();
}

/**
 * Delete centers.
 */
void IsoData::DeleteCenters()
{
	int c = centers.size(), i;

	for ( i = 0; i < c; i++ )
		delete centers[i];

	centers.clear();
}

/**
 * If points have been generated already, just set them to all points.
 */
void IsoData::SetPoints( const KMpointArray& all )
{
	for ( int i = 0; i < ImageSizeInByte; i++ )
	{
		allPoints[i] = new Point( NumBands, all[i] );
	}
}

/**********************************************************************/
/* Builds the 'filter; field based on the current points and centers. */
/**********************************************************************/
void IsoData::BuildKMfilterCenters()
{
	int i, s = samples.size();

	KMpointArray all = new KMpoint[s];

	for ( i = 0; i < s; i++ )
	{
		all[i] = kmAllocCopyPt( NumBands, allPoints[samples[i]]->getPoint() );
	}

	int c = centers.size();
	vector< KMpoint > cnts;

	for ( i = 0; i < c; i++ )
	{
		KMpoint P = new KMcoord[NumBands];
		kmCopyPt( NumBands, centers[i]->getPoint(), P );
		cnts.push_back( P );
	}

	data = new KMdata( NumBands, s, all );
	filter = new KMfilterCenters( c, *data, cnts );
	DeleteCenters();

}


//copies the image centers to filter centers.
void IsoData::SetFilterCenters()
{

	filter->SetFilterCenters( centers );
	DeleteCenters();

}

// copies the filter centers to the image centers.
void IsoData::SetImageCenters()
{

	vector< KMpoint > v = filter->getCtrPts();
	int num = v.size(), i;

	DeleteCenters();

	for ( i = 0; i < num; i++ )
	{
		Point* p = new Point();
		p = p->AllocPoint( NumBands, v[i] );
		centers.push_back( p );
	}

}

/**
 * Get total centers.
 */
int IsoData::GetNumCenters()
{
	return centers.size();
}

/******************************************************************************/
/* sampleCenters selects a set of initial centers from the set of allPoints.  */
/******************************************************************************/
void IsoData::SampleCenters()
{
	srand( sample_seed );

	int i, to_add;
	for ( i = 0; i < NumClusters; i++ )
	{
		to_add = rand() % ImageSizeInByte;

		//If the random point is not a duplicate, add it to the list of centers.
		Point* p;

		p = new Point( *allPoints[to_add] );

		if ( !p )
		{
			cout << "Error 7: Memory Allocation Failed" << endl;
			exit( 1 );
		}

		if ( FindCenter( p ) == -1 )
		{
			centers.push_back( p );
		}
		else
		{
			i--;
		}
	}
}

/**
 * Sample points randomly to perform iterative clustering on them.
 */
void IsoData::SamplePoints( double s )
{
	srand( sample_seed );

	int n = ImageSizeInByte;

	double probability = s / n;

	for ( int i = 0; ( i < ImageSizeInByte ) && ( s > 0 ); i++ )
	{
		//generate a random number between 0 and 1
		double x = ( (double) rand() / (double) ( RAND_MAX ) );

		if ( x <= probability )
		{
			samples.push_back( i );

			s = s - 1;
		}

		n = n - 1;
		probability = s / n;
	}
}

bool IsoData::WasDeleted()
{
	return Deleted;
}
/***********************************************************************************/
/* Given all the points in the image we take sampled points to determine	   */
/* cluster means, for having the final clasified image, one can run this function  */
/* with increment=1 to have all the points classified				   */
/***********************************************************************************/
void IsoData::PutInCluster()
{
	Deleted = false;

	if ( filter != NULL )
		filter->ComputeDist();
	else
		kmError( "In Image.cc :PutInCluster,KMfilterCenter object field is NULL.\nNeed to Build this object first", KMabort );

}
/*******************************************************************************/
/* searches in the 'centers' vector to see if it can find a particular point   */
/* It returns the index of the point in the vector if found, and -1 otherwise  */
/*******************************************************************************/
int IsoData::FindCenter( Point* to_find )
{
	int value = -1;
	int size = centers.size(), i;

	for ( i = 0; i < size; i++ )
	{
		if ( ( *centers[i] ) == ( *to_find ) )
		{
			value = i;
			break;
		}

	}

	return value;
}
/************************************************************************************/
/* This function checks number of points in each cluster. If number of points in    */
/* any cluster is less than NumSamples (desired minimum number of points in each    */
/* cluster), then that cluster and its center gets deleted (Note: only the cluster  */
/* from 'clusters' vector gets deleted, the points in the cluster do not get deleted*/
/* from set of sampled points for further iterative clustering.                     */
/************************************************************************************/
void IsoData::PostAnalyzeClusters( bool verbose )
{
	int i, num_clus, index, count;

	num_clus = clusters.size();

	vector< vector< int > >::iterator clusters_it = clusters.begin();

	count = 0;
	for ( i = 0, index = 0; index < num_clus; i++, clusters_it++, index++ )
	{
		if ( clusters[i].size() < NumSamples )
		{
			clusters.erase( clusters_it );
			Deleted = true;
			clusters_it--;
			i--;
			count++;
		} else
			centers_to_keep.push_back( index );
	}

	if ( Deleted )
	{
		if( verbose )
			std::cout << "\tDeleted " << count << " cluster(s)." << endl;

		filter->UpdateCenters( centers_to_keep );

	}

}

/****************************************************************************************/
/* Update each cluster center by setting it to the sample mean of its corresponding set */
/****************************************************************************************/
void IsoData::UpdateCenters()
{
	// moves each center to the centroid of the cluster
	filter->lloyd1Stage( centers_to_keep );
	centers_to_keep.clear();

}

/************************************************************************************/
/* we need to calculate the average of squared distances of all points in a cluster */
/* from that cluster center.  							    */
/***********************************************************************************/
void IsoData::CalculateAverageDistances()
{
	int num_centers = clusters.size(), i;

	average_distances = filter->getDists( true );

	for ( i = 0; i < num_centers; i++ )
	{
		// here we are taking square root of average of squared distances
		// in order to make this algorithm's result as close as possible to
		// the standard ISOCLUS.
		average_distances[i] = sqrt( average_distances[i] / clusters[i].size() );
	}
}

/************************************************************************************/
/* This function retunrs:                                                           */
/* over all average of squared distances of all points in the sample set form	    */
/* their closest center.		                                                    */
/************************************************************************************/
double IsoData::OverallAverageDistances()
{
	double sum = 0;
	int i, size = clusters.size();

	for ( i = 0; i < size; i++ )
		sum += ( clusters[i].size() ) * average_distances[i];

	OverallD = sum / samples.size();

	return OverallD;
}

/*****************************************************************************************/
/* This function calculates a vector for each cluster (a point for each cluster).        */
/* Each coordinate D of a particular  cluster's standard deviation vector is the         */
/* standard deviation of all coordinate D-s of all points in that cluster from coordinate*/
/* D of the same cluster's center.                                                       */
/*****************************************************************************************/
void IsoData::CalculateSTDVector()
{
	if ( filter != NULL )
	{
		STDVector = filter->getStdv( true );
	} else
	{
		kmError( "IsoData::CalculateSTDVector-Cannot get the standard deviation vector, filter is null", KMabort );
	}
}
/***************************************************************************/
/* calculates the maximum element in each column of STDVector. Since after */
/* this step we do not need values of STDVector, and each IsoClus iteration*/
/* needs to recalculate STDVector, we will delete the allocated memory for */
/* STDVector at this function as well.					   */
/**************************************************************************/
void IsoData::CalculateVmax()
{
	int c, b;
	int num_clusters = clusters.size();

	if ( Vmax != NULL )
	{
		delete[] Vmax;
		Vmax = NULL;
	}

	if ( Vmax_index != NULL )
	{
		delete[] Vmax_index;
		Vmax_index = NULL;
	}

	Vmax = new double[num_clusters];
	if ( !Vmax )
	{
		std::cerr << "Memory Allocation for 'Vmax' Failed." << endl;
		std::cerr << "Exitting the program..." << endl;
		exit( 1 );
	}

	Vmax_index = new int[num_clusters];
	if ( !Vmax_index )
	{
		std::cerr << "Memory Allocation for 'Vmax_index' Failed." << endl;
		std::cerr << "Exitting the program..." << endl;
		exit( 1 );
	}

	for ( c = 0; c < num_clusters; c++ )
	{
		Vmax[c] = STDVector[c][0];
		Vmax_index[c] = 0;

		for ( b = 1; b < NumBands; b++ )
		{
			if ( STDVector[c][b] > Vmax[c] )
			{
				Vmax[c] = STDVector[c][b];
				Vmax_index[c] = b;
			}
		}
	}
}
/**************************************************************************************/
/* This function evaluates 3 conditions X, A, and B and evaluates the value of:	      */
/* (X and (A or B)) for each cluster.						                          */
/* These conditions are described in step 10 of IsoClus algorithm.  See top of        */
/* IsoClus.cc for references.							                              */
/* It adds the cluster numbers of all those clusters that satisfied (X and (A or B))  */
/**************************************************************************************/
vector< int > IsoData::ShouldSplit( double stdv )
{
	vector< int > to_return;
	bool X = false, A = false, B = false;
	int num_clusters = clusters.size();

	//evaluateing X
	for ( int i = 0; i < num_clusters; i++ )
	{
		if ( Vmax[i] > stdv )
		{
			X = true;
		}
		//the goal is to evaluate [ X and (A or B) ]
		if ( X )
		{
			//evaluate A
			if ( ( average_distances[i] > OverallD ) && ( clusters[i].size() > 2 * ( NumSamples + 1 ) ) )
			{
				A = true;
			}
			if ( A )
				to_return.push_back( i );

			else
			{
				//evaluate B: compares the actual number of clusters with the desired number of
				//clusters.
				if ( num_clusters <= ( NumClusters / 2 ) )
				{

					B = true;
				}
				if ( B )
					to_return.push_back( i );

			}
		}

		// reset the booleans for next iteration
		X = false;
		A = false;
		B = false;
	}
	return to_return;
}

/****************************************************************************************/
/* It splits all those clusters whose index appears in 'to_split' parameter.  It does   */
/* so by calculating two new centers.  It assigns one of the centers to the same cluster*/
/* that is being split, and add the other new center to the end of 'centers' array.     */
/****************************************************************************************/
void IsoData::Split( vector< int > to_split, bool verbose )
{
	int i, size = to_split.size(), j, index, num_centers = centers.size();
	double Gj, current;
	Point Zminus, Zplus;
	try
	{
		for ( j = 0; j < size && num_centers < std::numeric_limits< int >::max(); j++, num_centers++ )
		{
			index = to_split[j];
			Gj = K*Vmax[index];

			Zminus = ( *centers[index] );
			Zplus = ( *centers[index] );

			//i would be the coordinate of STDV vector which had the maximum component.

			i = Vmax_index[index];
			current = Zminus.getCoordinate( i );

			Zminus.setCoordinate( i, current - Gj );
			Zplus.setCoordinate( i, current + Gj );

			// this erases the previous center value and updates it by what is called
			// Zplus in in ISOCLUS algorithm step 10, see references in IsoClus.cc

			( *centers[index] ) = Zminus;
			centers.push_back( new Point( Zplus ) );

			if ( verbose )
				std::cout << "\tSplit cluster " << index + 1 << "." << endl;

		}
	}
	catch ( int e )
	{
		std::cerr << "Exception occured in IsoData::Split function: " << e << endl;
		exit( 1 );
	}
}

/**
 * Compute center distances.
 */
void IsoData::ComputeCenterDistances()
{
	unsigned int i, j, size, count;
	double dist;
	bool Emptied = false;

	size = centers.size();

	if ( size != CenterDistances.size() )
	{
		CenterDistances.clear();
		Emptied = true;
	}

	if ( size == 0 || Emptied )
	{
		for ( i = 0; i < size - 1; i++ )
			for ( j = i + 1; j < size; j++ )
			{
				dist = sqrt( centers[i]->Norm2DistanceSquared( centers[j] ) );
				PairDistanceNode n =
				{ dist, i, j };
				CenterDistances.push_back( n );

			}
	}

	else //just need to update "dist" values rather than allocating all structure nodes
	{
		//note that 'count' is suppose to be the index number of the next element in
		//CenterDistances array.
		for ( i = 0, count = 0; i < size - 1; i++ )
		{
			for ( j = i + 1; j < size; j++, count++ )
			{
				dist = sqrt( centers[i]->Norm2DistanceSquared( centers[j] ) );
				CenterDistances[count].dist = dist;
			}
		}
	}
}

/************************************************************************************/
/* Searches the list of pair distance nodes and selects those pairs whose distances */
/* from each other are less than parameter LUMP.  Orders the list in ascending order*/
/* and selects the first MAXPAIR pairs as candid pairs for lumping.		    */
/************************************************************************************/
vector< PairDistanceNode > IsoData::FindLumpCandidates( double lump, int MAXPAIR )
{
	int count = 0;
	int size = CenterDistances.size();
	vector< PairDistanceNode > lump_candidates;

	// sort the list of center pairs based on '<' operator which is overloaded to compare
	// 'dist' field of each node. (sort based on distances)
	sort( CenterDistances.begin(), CenterDistances.end() );

	for ( int i = 0; i < size && lump != 0 && count < MAXPAIR; i++ )
	{
		if ( CenterDistances[i] < lump )
		{

			lump_candidates.push_back( CenterDistances[i] );
			count += 1;
		}
	}

	//If there are more candidates to lump, than the MAXPAIR, select the first MAXPAIR candidates

	if ( count > MAXPAIR )
		lump_candidates.erase( lump_candidates.begin() + MAXPAIR, lump_candidates.end() );

	return lump_candidates;

}
/**********************************************************************************/
/* Eligible clusters among to_lump vector will be lumped.                         */
/* Each cluster/center can be lumped only once, and thus not all the centers      */
/* associated with elements of to_lump vector will be lumped.                     */
/**********************************************************************************/
void IsoData::Lump( const vector< PairDistanceNode >& to_lump, bool verbose )
{
	int to_lump_size = to_lump.size();
	int orig_centers_size = centers.size();
	int i, clus1_size, clus2_size, clus1, clus2, used_index, count = 0;

	try
	{
		bool* used_centers = new bool[orig_centers_size];

		for ( i = 0; i < orig_centers_size; i++ )
			used_centers[i] = 0;

		for ( i = 0; i < to_lump_size; i++ )
		{
			clus1 = to_lump[i].c1;
			clus2 = to_lump[i].c2;

			// only go about lumping the two clusters if neither of the centers had been used before.
			if ( !used_centers[clus1] && !used_centers[clus2] && ( orig_centers_size - count ) > MINCLUS )
			{
				//calculate the new center
				clus1_size = clusters[clus1].size();
				clus2_size = clusters[clus2].size();

				Point new_center = ( ( *centers[clus1] ) * clus1_size ) + ( ( *centers[clus2] ) * clus2_size );
				new_center = new_center / ( clus1_size + clus2_size );
				centers.push_back( new Point( new_center ) );

				//merge the two clusters into a new one, and add the new cluster to the end of vector of clusters.
				vector< int > new_cluster;
				new_cluster.insert( new_cluster.end(), clusters[clus1].begin(), clusters[clus1].end() );
				new_cluster.insert( new_cluster.end(), clusters[clus2].begin(), clusters[clus2].end() );
				clusters.push_back( new_cluster );

				used_centers[clus1] = 1;
				used_centers[clus2] = 1;

				if( verbose )
					std::cout << "\tLumped clusters " << clus1 + 1 << " and " << clus2 + 1 << "." << endl;

				count++;
			} // if could be lumped
		} // for loop

		vector< Point* >::iterator centers_it = centers.begin();
		vector< vector< int > >::iterator clusters_it = clusters.begin();

		for ( i = 0, used_index = 0; used_index < orig_centers_size; i++, centers_it++, clusters_it++, used_index++ )
		{
			//if this cluster has been lumped remove it and its center.
			if ( used_centers[used_index] )
			{
				clusters.erase( clusters_it );
				delete centers[i];
				centers.erase( centers_it );
				i--;
				centers_it--;
				clusters_it--;

			} // if
		} // for
	} //try
	catch ( int e )
	{
		std::cerr << "Exception occured in IsoData::Lump function: " << e << endl;
		exit( 1 );
	}
}

/***************************************************************************************/
/* Select all points for the last iterative clustering.                                */
/***************************************************************************************/
void IsoData::PreFinalClustering()
{
	//since it is going to be the last round, we need to classify every single point

	samples.clear();
	for ( int i = 0; i < ImageSizeInByte; i++ )
	{
		samples.push_back( i );
		if ( allPoints[i] == NULL )
		{
			std::cerr << "*** ERROR ***: one of the points: NULL!" << std::endl;
		}
	}
	// since it is the last iteration, build a new KMfilterCenters to consider all points
	// and the latest centers.
	if ( filter != NULL )
	{
		delete filter;
		filter = NULL;
	}
	if ( data != NULL )
	{
		delete data;
		data = NULL;
	}
	BuildKMfilterCenters();
}

/****************************************************************************************/
/* Print a report of current clustering results. (Number of clusters, their, sizes, etc */
/****************************************************************************************/
void IsoData::GenerateReport( ostream* out )
{

	int i, num_clusters = clusters.size();
	*out << "\n\tNumber of Clusters: " << clusters.size() << endl;
	*out << "\t===============================================================================\n";
	for ( i = 0; i < num_clusters; i++ )
	{
		*out << "\tCluster " << i + 1 << ":\n\tSize: " << clusters[i].size();
		*out << "\n\tAverage of squared distances: ";
		*out << average_distances[i] << endl;
		*out << "\tCenter:  ";
		centers[i]->print( out );
		*out << "\t===============================================================================\n";

	}
	*out << "\tOverall average of squared distances of all points from their cluster center: " << OverallD << endl;
}


