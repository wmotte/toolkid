#include "Point.h"
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <numeric>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNormalizeImageFilter.h"
#include "itkMaskImageFilter.h"

// typedefs...
typedef itk::Image< float, 3 > ImageType;
typedef itk::Image< unsigned char, 3 > MaskType;
typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::NormalizeImageFilter< ImageType, ImageType > NormalizeImageFilterType;
typedef itk::MaskImageFilter< ImageType, MaskType, ImageType > MaskImageFilterType;

class IsoDataImage {

	// ***********************************************************************************
protected:

	std::vector< Point > points;
	std::vector< Point > clusters;

	std::vector< std::vector< double > > distances; // rows = points; columns   = Points2ClusterDistance
	std::vector< std::vector< Point > > clusterMembers; // rows = clusters; columns = PointsInCluster.
	std::vector< std::vector< double > > stdVectors;

	std::vector< unsigned int > vmax_index;
	std::vector< double > vmax;
	std::vector< double > averageDistances; // step 5...

	double overallAverageDistance; // step 6...
	unsigned int kinit;
	unsigned int dimensions;

	/**
	 * Return random sample from points.
	 */
	unsigned int getRandomSampleNumber() {

		if ( points.empty() ) {
			std::cerr << "Points are not initialized!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// get random number in range of points...
		unsigned int random = (unsigned) ( rand() % ( points.size() - 1 ) );

		// sanity check...
		if ( random < points.size() ) {
			return random;
		} else {
			return getRandomSampleNumber();
		}
	}

	/**
	 * Read input and convert to 1D float channel.
	 */
	void inputToChannel( const std::string& inputFileName, std::vector< float >& channel, const std::string& mask, bool normalize ) {

		// read input...
		ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
		reader -> SetFileName( inputFileName );

		ImageType::Pointer input;

		try {
			if ( normalize ) {
				NormalizeImageFilterType::Pointer filter = NormalizeImageFilterType::New();
				filter -> SetInput( reader -> GetOutput() );
				input = filter -> GetOutput();
				filter -> Update();
			} else {
				input = reader -> GetOutput();
				reader -> Update();
			}
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "*** Error *** Could not read: " << input << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}

		ImageType::RegionType region = input -> GetLargestPossibleRegion();
		IteratorType it( input, region );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) {

			ImageType::IndexType index = it.GetIndex();
			ImageType::PixelType pixel = input -> GetPixel( index );

			channel.push_back( pixel );
		}
	}

	/**
	 * Convert channels to points array.
	 */
	void channelsToPoints( const std::vector< std::vector< float > >& channels, std::vector< Point >& points ) {

		unsigned int dimensions = channels.size();
		unsigned int num_voxels = channels[0].size();

		// for each voxel...
		for ( unsigned int i = 0; i < num_voxels; i++ ) {

			std::vector< double > voxels;

			// extract all values...
			for ( unsigned int j = 0; j < dimensions; j++ ) {

				double value = channels[j][i];

				voxels.push_back( value );
			}

			Point point = Point( voxels );
			points.push_back( point );
		}
	}

	/**
	 * Calculate the distances of all points from centers, and look up these distances later.
	 */
	void calculateDistances() {

		distances.clear();

		for ( unsigned int i = 0; i < points.size(); i++ ) {

			std::vector< double > distance2clusters;

			for ( unsigned int j = 0; j < clusters.size(); j++ ) {

				double value = points[i].norm2DistanceSquared( clusters[j] );
				distance2clusters.push_back( value );

			}
			distances.push_back( distance2clusters );
		}
	}

	/**
	 * Return the lowest element: O(n).
	 */
	double findLowestValue( const std::vector< double >& values ) {
		std::vector< double >::const_iterator lowest = std::min_element( values.begin(), values.end() );
		return *lowest;
	}

	/**
	 * Return the lowest element index: O(n).
	 */
	unsigned int findLowestValueIndex( const std::vector< double >& values ) {
		std::vector< double >::const_iterator lowest = std::min_element( values.begin(), values.end() );
		return ( lowest - values.begin() );
	}

	/**
	 * Set value to highest element: O(n). Return index.
	 */
	int findHighestValue( const std::vector< double >& values, double& value ) {
		std::vector< double >::const_iterator highest = std::max_element( values.begin(), values.end() );
		value = *highest;
		return ( highest - values.begin() );
	}

	/**
	 * Return average point.
	 */
	Point averagePoint( std::vector< Point >& points ) {

		unsigned int size = points.size();

		// sanity check...
		if ( size == 0 ) {
			std::cerr << "Points are empty! Not able to construct average point!" << std::endl;
			exit( EXIT_FAILURE );
		}

		Point average = points[0];

		for ( unsigned int i = 1; i < points.size(); i++ ) {
			average = averagePoint( average, points[i] );
		}

		return ( dividePoint( average, size ) );
	}

	/**
	 * Return average Point.
	 */
	Point averagePoint( Point& one, Point& two ) {

		if ( one.getDimension() != two.getDimension() ) {
			std::cerr << "Dimensions are not similar in averagePoint function!" << std::endl;
			exit( EXIT_FAILURE );
		}

		std::vector< double > values1 = one.getValues();
		std::vector< double > values2 = two.getValues();
		std::vector< double > final( values1.size() );

		for ( unsigned int i = 0; i < values1.size(); i++ ) {
			final[i] = ( values1[i] + values2[i] ) / 2.0;
		}

		return Point( final );
	}

	/**
	 * Divide all values in Point with constant.
	 */
	Point dividePoint( Point point, float constant ) {

		if ( constant == 0.0 ) {
			std::cerr << "Cannot divide by 0!" << std::endl;
			exit( EXIT_FAILURE );
		}

		std::vector< double > values = point.getValues();
		std::vector< double > final( values.size() );

		for ( unsigned int i = 0; i < values.size(); i++ ) {
			final[i] = ( values[i] / constant );
		}

		return Point( final );
	}

	// ***********************************************************************************
public:
	static const double K = 0.5;
	static const unsigned int MINCLUS = 1;

	struct PairDistanceNode {
		double dist;
		unsigned int c1;
		unsigned int c2;
		bool operator<( PairDistanceNode to_compare ) const {
			return ( dist < to_compare.dist );
		}
		bool operator<( double to_compare ) const {
			return ( dist < to_compare );
		}
	};

	/**
	 * Main constructor. Responsible for loading images in memory.
	 */
	IsoDataImage( const std::vector< std::string >& inputs, const std::string& output, const std::string& mask, bool normalize ) {

		// set dimensions...
		dimensions = inputs.size();

		// set seed when class loads... (for random cluster selection).
		srand( (unsigned) time( 0 ) );

		std::vector< std::vector< float > > channels;

		// convert image list to floats...
		for ( unsigned int i = 0; i < inputs.size(); i++ ) {
			std::vector< float > channel;
			inputToChannel( inputs[i], channel, mask, normalize );
			channels.push_back( channel );
		}

		channelsToPoints( channels, points );
	}

	/**
	 * Initialize clusters with random points taken from data.
	 */
	void initClusters( unsigned int num_clusters ) {

		kinit = num_clusters;

		clusters.clear();

		for ( unsigned int i = 0; i < num_clusters; i++ ) {
			Point randomPoint = points[getRandomSampleNumber()];
			clusters.push_back( randomPoint );
		}
	}

	/**
	 * With calculated distances we are able to assign each point
	 * to the closest cluster.
	 */
	void assignPointsToClosestClusterCenter() {

		// first calculate distances from points to cluster centers...
		calculateDistances();

		// initialize multi-array output...
		std::vector< std::vector< Point > > output( clusters.size() );

		// for each point...
		for ( unsigned int i = 0; i < points.size(); i++ ) {

			// get distances to all cluster centers...
			std::vector< double > distance = distances[i];

			// calculate index for closest cluster center...
			unsigned int index = findLowestValueIndex( distance );

			// assign point to closest cluster center...
			output[index].push_back( points[i] );
		}

		clusterMembers = output;
	}

	/**
	 * Discard clusters with fewer than minimal required points per cluster...
	 * If there is a minimal clusters, than return true.
	 */
	bool discardSmallClusters( unsigned int minPointsPerCluster ) {

		bool clusterRemoved = false;

		// new clusterMembers collection...
		std::vector< std::vector< Point > > output;

		for ( std::vector< std::vector< Point > >::iterator it = clusterMembers.begin(); it != clusterMembers.end(); ++it ) {
			unsigned int size = ( *it ).size();

			if ( size > minPointsPerCluster ) {
				output.push_back( *it );
			} else {
				clusterRemoved = true; // one or more removed clusters: so return value is true.
			}
		}

		clusterMembers = output;

		return clusterRemoved;
	}

	/**
	 * Update each cluster center: sum all points in cluster and divide by total points in cluster...
	 */
	void updateClusters() {

		clusters.clear();

		for ( unsigned int i = 0; i < clusterMembers.size(); i++ ) {
			Point cluster = averagePoint( clusterMembers[i] );
			clusters.push_back( cluster );
		}
	}

	/**
	 * Compute average distance of samples in clusters.
	 */
	void computeAverageDistance() {

		// get latest distances to cluster centers...
		calculateDistances();

		// clean previous values...
		averageDistances.clear();

		for ( unsigned int i = 0; i < clusters.size(); i++ ) {
			std::vector< double > distance = distances[i];

			double sum = std::accumulate( distance.begin(), distance.end(), 0 );
			double size = distance.size();

			averageDistances.push_back( sum / size );
		}
	}

	/**
	 * Compute overall average distance of samples to cluster centers.
	 */
	void computeOverallAverageDistance() {

		// check empty distances
		if ( averageDistances.empty() ) {
			std::cerr << "Distances are empty! Not able to construct overall average distance!" << std::endl;
			exit( EXIT_FAILURE );
		}

		double sum = std::accumulate( averageDistances.begin(), averageDistances.end(), 0 );
		double size = averageDistances.size();

		overallAverageDistance = ( sum / size );
	}

	/**
	 * Return number of cluster centers...
	 */
	unsigned int getNumCenters() {
		return clusters.size();
	}

	void calculateSTDVector() {

		stdVectors.clear();

		// for each cluster...
		for ( unsigned int cluster_index = 0; cluster_index < clusters.size(); cluster_index++ ) {

			std::vector< Point > clusterPoints = clusterMembers[cluster_index];
			unsigned int cluster_size = clusterPoints.size();

			std::vector< double > stdVector;

			// for each dimension...
			for ( unsigned int dim = 0; dim < dimensions; dim++ ) {

				double total = 0.0;
				double avg = 0.0;

				// for each point...
				for ( unsigned int point_index = 0; point_index < cluster_size; point_index++ ) {

					Point point = clusterPoints[point_index];
					Point cluster = clusters[cluster_index];

					total += pow( ( point.getCoordinate( dim ) - cluster.getCoordinate( dim ) ), 2 );
				}

				if ( cluster_size != 0 ) {
					avg = ( total / cluster_size );
				}

				stdVector.push_back( sqrt( avg ) );
			}

			stdVectors.push_back( stdVector );
		}
	}

	/**
	 * Find the maximum element of each column of stdVectors.
	 */
	void calculateVmax() {

		vmax.clear();
		for ( unsigned int i = 0; i < stdVectors.size(); i++ ) {
			double value;
			int index = findHighestValue( stdVectors[i], value );
			vmax.push_back( value );
			vmax_index.push_back( index );
		}
	}

	/**
	 * The vector to_split will contain integers that represent the cluster numbers that need to be split.
	 */
	std::vector< unsigned int > shouldSplit( double STDV, unsigned int minPointsPerCluster ) {

		std::vector< unsigned int > to_return;
		bool X = false;
		bool A = false;
		bool B = false;

		// Evaluating X...
		for ( unsigned int i = 0; i < clusters.size(); i++ ) {
			if ( vmax[i] > STDV ) {
				X = true;
			}

			// Evaluate [ X and (A or B) ] ...
			if ( X ) {
				// evaluate A...
				if ( ( averageDistances[i] > overallAverageDistance ) && ( clusterMembers[i].size() > 2 * ( minPointsPerCluster + 1 ) ) ) {
					A = true;
				}

				if ( A ) {
					to_return.push_back( i );
				}

				else {
					//evaluate B: compares the actual number of clusters with the desired number of
					//clusters.
					if ( clusters.size() <= ( kinit / 2 ) ) {
						B = true;
					}

					if ( B ) {
						to_return.push_back( i );
					}
				}
			}

			// reset the booleans for next iteration
			X = false;
			A = false;
			B = false;
		}

		return to_return;
	}

	/**
	 * Split all those clusters whose index appears in 'to_split' parameter.
	 * It does so by calculating two new centers.
	 * It assigns one of the centers to the same cluster that is being split,
	 * and add the other new center to the end of 'centers' array.
	 */
	void split( const std::vector< unsigned int >& to_split ) {

		for ( unsigned int i = 0; i < to_split.size(); i++ ) {

			double Gj = K * vmax[i];

			Point Zminus = clusters[i];
			Point Zplus = clusters[i];

			// Get coordinate of STDV vector which had the maximum component...
			double current = Zminus.getCoordinate( vmax_index[i] );

			Zminus.setCoordinate( i, current - Gj );
			Zplus.setCoordinate( i, current + Gj );

			// replace old cluster and add additional one...
			clusters[i] = Zminus;
			clusters.push_back( Zplus );
		}
	}

	/**
	 * Computes the pairwise distance between cluster centers.
	 */
	std::vector< PairDistanceNode > computeCenterDistances() {

		std::vector< PairDistanceNode > centerDistances;

		for ( unsigned int i = 0; i < clusters.size() - 1; i++ ) {
			for ( unsigned int j = i + 1; j < clusters.size(); j++ ) {
				double dist = sqrt( clusters[i].norm2DistanceSquared( clusters[j] ) );
				PairDistanceNode n = { dist, i, j };
				centerDistances.push_back( n );
			}
		}

		return centerDistances;
	}

	/**
	 * Searches the list of pair distance nodes and selects those pairs whose distances
	 * from each other are less than parameter LUMP.  Orders the list in ascending order
	 * and selects the first MAXPAIR pairs as candid pairs for lumping.
	 */
	std::vector< PairDistanceNode > findLumpCandidates( std::vector< PairDistanceNode >& centerDistances, double LUMP, unsigned int MAXPAIR ) const {

		unsigned int count = 0;
		unsigned int size = centerDistances.size();
		std::vector< PairDistanceNode > lump_candidates;

		// sort the list of center pairs based on '<' operator which is overloaded to compare
		// 'dist' field of each node. (sort based on distances)...
		std::sort( centerDistances.begin(), centerDistances.end() );

		for ( unsigned int i = 0; i < size && ( LUMP != 0 ) && ( count < MAXPAIR ); i++ ) {
			if ( centerDistances[i] < LUMP ) {
				lump_candidates.push_back( centerDistances[i] );
				count += 1;
			}
		}

		sort( lump_candidates.begin(), lump_candidates.end() );

		// if there are more candidates to lump, than the MAXPAIR,
		// select the first MAXPAIR candidates...
		if ( count > MAXPAIR )
			lump_candidates.erase( lump_candidates.begin() + MAXPAIR, lump_candidates.end() );

		return lump_candidates;
	}

	/**
	 * Eligible clusters among to_lump vector will be lumped.
	 * Each cluster/center can be lumped only once, and thus not all the centers
	 * associated with elements of to_lump vector will be lumped.
	 */
	void lump( std::vector< PairDistanceNode > to_lump ) {

		unsigned int to_lump_size = to_lump.size();
		unsigned int orig_centers_size = clusters.size();
		unsigned int clus1_size, clus2_size, clus1, clus2, used_index, count = 0;

		bool* used_centers = new bool[orig_centers_size];

		for ( unsigned int j = 0; j < orig_centers_size; j++ ) {
			used_centers[j] = 0;
		}

		for ( unsigned int k = 0; k < to_lump_size; k++ ) {
			clus1 = to_lump[k].c1;
			clus2 = to_lump[k].c2;

			// only go about lumping the two clusters if neither of the centers had been used before.
			if ( !used_centers[clus1] && !used_centers[clus2] && ( orig_centers_size - count ) > MINCLUS ) {
				//calculate the new center
				clus1_size = clusterMembers[clus1].size();
				clus2_size = clusterMembers[clus2].size();

				Point new_center = ( clusters[clus1] * clus1_size ) + ( ( clusters[clus2] ) * clus2_size );
				new_center = new_center / ( clus1_size + clus2_size );
				clusters.push_back( new_center );

				used_centers[clus1] = 1;
				used_centers[clus2] = 1;
				count++;
			}
		}

		delete[] used_centers;

	}

	void writeOutput( const std::string& output ) {
		// TODO
	}

	// ****************************************************************************
	// DEBUG functions...
	// ****************************************************************************

	/**
	 * Print cluster values (debug).
	 */
	void printClusters() {
		for ( unsigned int i = 0; i < clusters.size(); i++ ) {
			std::cout << "Value: " << i << " -> " << clusters[i].getStringValues() << std::endl;
		}
	}

	/**
	 * Print cluster sizes (debug).
	 */
	void printClusterSizes() {
		for ( unsigned int i = 0; i < clusterMembers.size(); i++ ) {
			std::cout << "Total members in cluster: " << i << " -> " << clusterMembers[i].size() << std::endl;
		}
	}

	/**
	 * Print cluster distances for given point index (debug).
	 */
	void printDistances( unsigned int index ) {
		std::cout << "Printing cluster distances for point index: " << index << std::endl;
		if ( index < distances.size() ) {
			std::vector< double > values = distances[index];
			std::cout << "Total cluster distances for given point are: " << values.size() << std::endl;

			for ( unsigned int i = 0; i < values.size(); i++ ) {
				std::cout << "Cluster (" << i << ") distance: " << values[i] << std::endl;
			}

			double lowestValue = findLowestValue( values );
			int lowestValueIndex = findLowestValueIndex( values );
			std::cout << "Lowest value is: " << lowestValue << std::endl;
			std::cout << "Lowest value index is: " << lowestValueIndex << std::endl;

		} else {
			std::cerr << "Index out of range" << std::endl;
		}
	}

};

