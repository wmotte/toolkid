#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkTimeProbe.h"
#include "IsoDataImage.cpp"

/**
 * ISODATA.
 */
class IsoData {

protected:

	// Split function; return true if we should split...
	bool split( IsoDataImage& img, double STDV, unsigned int minPointsPerCluster ) {

		bool split = false;

		//STEP 7:
		img.calculateSTDVector();

		//STEP 8:
		img.calculateVmax();

		// STEP 9:
		// the vector to_split will contain integers that
		// represent the cluster numbers that need to be split.

		std::vector< unsigned int > to_split = img.shouldSplit( STDV, minPointsPerCluster );

		if ( to_split.size() != 0 ) {
			img.split( to_split );
			split = true;
		}

		return split;
	}

	// Merge...
	void merge( IsoDataImage& img, double LUMP, int MAXPAIR ) {
		// STEP 11:
		std::vector< IsoDataImage::PairDistanceNode > centerDistances = img.computeCenterDistances();

		// STEP 12:
		std::vector< IsoDataImage::PairDistanceNode > to_lump = img.findLumpCandidates( centerDistances, LUMP, MAXPAIR );

		// STEP 13:
		if ( to_lump.size() != 0 ) {
			img.lump( to_lump );
		}
	}

	// ************************************************************************************************
public:

	/**
	 * Run.
	 */
	void run( const std::vector< std::string >& inputs, const std::string& output,
			const std::string& mask, unsigned int NUMCLUS, unsigned int SAMPRM,
			unsigned int MAXITER, double STDV, double LUMP, int MAXPAIR, bool normalize ) {

		itk::TimeProbe probe;
		probe.Start();

		// initialize raw data...
		IsoDataImage img = IsoDataImage( inputs, output, mask, normalize );

		// STEP 1: arbitrarily choose k and init clusters...
		img.initClusters( NUMCLUS );

		for ( unsigned int i = 0; i < MAXITER; i++ ) {

			std::cout << "Iter: "  << i << std::endl;
			// STEP 2: assign each point to the clostest cluster center...
			img.assignPointsToClosestClusterCenter();

			// STEP 3: discard clusters with less than min. samples per cluster...
			// If any clusters were deleted, then go back to STEP 2...
			if ( img.discardSmallClusters( SAMPRM ) ) {
				continue;
			}

			// STEP 4: update each remaining cluster center...
			img.updateClusters();

			// STEP 5: compute the average distance of points in clusters from
			// their corresponding cluster center...
			// also compute overall average distance...
			img.computeAverageDistance();
			img.computeOverallAverageDistance();

			// STEP 6:
			if ( i == MAXITER ) {
				LUMP = 0.0;
				// goto STEP 9:
				merge( img, LUMP, MAXPAIR );
			} else if ( ( i % 2 == 0 ) || ( img.getNumCenters() >= 2 * NUMCLUS ) ) {
				// goto STEP 9:
				merge( img, LUMP, MAXPAIR );
			} else if ( img.getNumCenters() <= ( NUMCLUS / 2 ) ) {
				// goto STEP 7:
				// if we are in last iteration and split is performed, we need to rerun from STEP 2 again
				// as cluster centers need to be updated before final run...
				if ( ( split( img, STDV, SAMPRM ) ) && ( i == MAXITER ) ) {
					i--;
				}
			}
		}

		img.writeOutput( output );

		probe.Stop();
		std::cout << "Total runtime: " << probe.GetMeanTime() << " s." << std::endl;
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::vector< std::string > inputs;
	std::string output;
	std::string mask;

	int numberOfClusters = 10; // kinit
	int minNumberOfClusterPoints = 30; // nmin			1/5 average kluster size -> nmin = n/5 kinit.
	int maxNumberOfIter = 20; // Imax			default 20
	double maxStdev = 0.1; // STDV			default 2 * sigma
	double minRequiredDistance = 0.001; // Lmin			default 0.001
	int maxPair = 4; // MAXPAIR
	bool normalize = false;	// normalize input images...

	tkd::CmdParser parser( argv[0], "Isodata" );

	parser.AddArgument( inputs, "inputs" ) -> AddAlias( "i" ) -> SetInput( "<strings>" ) -> SetDescription( "Input images" ) -> SetRequired(
			true ) -> SetMinMax( 1, 10000 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output cluster image file" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( mask, "mask" ) -> AddAlias( "m" ) -> SetInput( "<string>" ) -> SetDescription( "Mask image" ) -> SetRequired(
				false ) -> SetMinMax( 0, 1 );

	parser.AddArgument( numberOfClusters, "clusters" ) -> AddAlias( "c" ) -> SetInput( "<int>" ) -> SetDescription(
			"Initial number of clusters (NUMCLUS)" ) -> SetMinMax( 1, 1 );

	parser.AddArgument( minNumberOfClusterPoints, "min-clusters" ) -> AddAlias( "m" ) -> SetInput( "<int>" ) -> SetDescription(
			"Minimum number of points that can form a cluster (SAMPRM)" ) -> SetMinMax( 1, 1 );

	parser.AddArgument( maxNumberOfIter, "max-itererations" ) -> AddAlias( "it" ) -> SetInput( "<int>" ) -> SetDescription(
			"Maximum number of iterations (MAXITER)" ) -> SetMinMax( 1, 1 );

	parser.AddArgument( maxStdev, "standard-deviation" ) -> AddAlias( "s" ) -> SetInput( "<double>" ) -> SetDescription(
			"Maximum standard deviation of points from their cluster center along each axis (STDV)" ) -> SetMinMax( 1, 1 );

	parser.AddArgument( minRequiredDistance, "min-distance" ) -> AddAlias( "d" ) -> SetInput( "<double>" ) -> SetDescription(
			" Minimum required distance between two cluster centers (LUMP)" ) -> SetMinMax( 1, 1 );

	parser.AddArgument( maxPair, "max-pairs" ) -> AddAlias( "p" ) -> SetInput( "<int>" ) -> SetDescription(
			"Maximum number of cluster pairs that can be merged per iteration (MAXPAIR)" ) -> SetMinMax( 1, 1 );

	parser.AddArgument( normalize, "normalize" ) -> AddAlias( "n" ) -> SetInput( "<bool>" ) -> SetDescription(
				"Normalize input images before clustering (default: false)" );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	IsoData isoData = IsoData();

	isoData.run( inputs, output, mask, numberOfClusters, minNumberOfClusterPoints, maxNumberOfIter, maxStdev, minRequiredDistance, maxPair, normalize );
}

