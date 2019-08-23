#include "tkdCmdParser.h"
#include "annCommon.h"
#include "graphCommon.h"
#include "annRelabel.h"

//Structure of the command line parameters
struct parameters
{
	std::string inputFileName;
	std::string outputFileName;
	std::string maskFileName;
	double threshold;
	int clusters;
	int repeat;
	bool verbose;
	bool noThreshold;
};

/**
 * Spectral clustering.
 */
class AnnClustering
{
public:

	typedef double PixelType;
	typedef itk::Image< unsigned short, 3 > LabelImageType;
	//typedef itk::Image< unsigned int, 3 > LabelImageType;

	/**
	 * Run.
	 */
	void Run( const parameters& list )
	{
		// if input is 2D image...
		if ( !ann::Ann< PixelType >::GetImageDimensions( list.inputFileName ) == 2 )
		{
			std::cerr << "*** ERROR ***: Dimensions of input image is not equal to 2 (matrix)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		graph::Graph< PixelType >::ImagePointerType input;
		graph::Graph< PixelType >::ImagePointerType output;

		graph::Graph< PixelType >::ReadMatrix( input, list.inputFileName );

		// threshold input matrix ...
		if( list.noThreshold )
			graph::Graph< PixelType >::GetMatrix( input, output, 0.0, std::numeric_limits< double >::min(), false, true );
		else
		{
			if ( list.verbose )
				std::cout << "Threshold matrix using fixed K..." << std::endl;
			graph::Graph< PixelType >::GetMatrix( input, output, 0.0, list.threshold, true, true );
		}

		// calculate laplacian ...
		if ( list.verbose )
			std::cout << "Calculate laplacian..." << std::endl;

		graph::Graph< PixelType >::MatrixType G;
		graph::Graph< PixelType >::GetMatrix( G, output );
		graph::Graph< PixelType >::MatrixType L;
		graph::Graph< PixelType >::DiagMatrixType D;
		graph::Graph< PixelType >::CalculateLaplacian( L, G, D, false );

		// solve eigen system ...
		if ( list.verbose )
			std::cout << "Solve generalized eigensystem..." << std::endl;

		graph::Graph< PixelType >::MatrixType eigenVectors;
		graph::Graph< PixelType >::VectorType eigenValues;
		graph::Graph< PixelType >::GeneralizedEigenCalculation( L, D, eigenVectors, eigenValues );

		if ( list.verbose )
			std::cout << "Kmeans++ clustering..." << std::endl;

		Cluster( G, eigenVectors, list );
	}

protected:

	/**
	 * Cluster.
	 */
	void Cluster( const graph::Graph< PixelType >::MatrixType& G, const graph::Graph< PixelType >::MatrixType& eigenVectors,
			const parameters& list )
	{

		std::vector< vnl_vector< PixelType > > totalClusterings;

		for ( int i = 0; i < list.repeat; i++ )
		{
			graph::Graph< PixelType >::VectorType clusters;
			graph::Graph< PixelType >::KMeansPP( eigenVectors, list.clusters, clusters );
			totalClusterings.push_back( clusters );
		}

		RelabelAndWrite( totalClusterings, list.clusters, list.maskFileName, list.outputFileName );
	}

	/**
	 * Relabel
	 */
	void RelabelAndWrite( const std::vector< vnl_vector< PixelType > >& matrix,
			unsigned int nLabels, const std::string& maskFileName,
			const std::string& outputFileNameRelabel )
	{
		if ( matrix.empty() )
		{
			std::cerr << "*** ERROR ***: could not calculate entropy (matrix empty)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		graph::Graph< PixelType >::Image4DPointerType image4D;
		graph::Graph< PixelType >::Zero4DImage( image4D, maskFileName, matrix.size() );
		graph::Graph< PixelType >::Image3DPointerType ref3D;

		// construct ref image...
		graph::Graph< PixelType >::FillImage( ref3D, maskFileName, matrix[0] );
		graph::Graph< PixelType >::Insert3DImageIn4DImage( image4D, ref3D, 0 );
		Relabel relabel;
		LabelImageType::Pointer refLabel3D;
		relabel.ConvertImageToLabelImage( refLabel3D, ref3D );

		for ( unsigned int i = 1; i < matrix.size(); i++ ) // skip first (=ref)...
		{
			graph::Graph< PixelType >::Image3DPointerType image3D;
			graph::Graph< PixelType >::FillImage( image3D, maskFileName, matrix[i] );

			LabelImageType::Pointer label3D;
			relabel.ConvertImageToLabelImage( label3D, image3D );
			LabelImageType::Pointer relabeled3D;
			relabel.RelabelImage( relabeled3D, label3D, refLabel3D );

			relabel.ConvertLabelImageToImage( image3D, relabeled3D );

			graph::Graph< PixelType >::Insert3DImageIn4DImage( image4D, image3D, i );
		}

		// Write relabeled file to output for check...
		graph::Graph< PixelType >::Write4DImage( image4D, outputFileNameRelabel );
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Cluster matrix using spectral clustering" );

	parameters list;

	list.clusters = 5;
	list.threshold = 50;

	list.repeat = 10;
	list.verbose = false;
	list.noThreshold = true;

	p.AddArgument( list.inputFileName, "input" )->AddAlias( "i" )->SetDescription( "Input file name (2D matrix)" )->SetRequired( true )->SetMinMax(
			1, 1 );
	p.AddArgument( list.maskFileName, "mask" )->AddAlias( "m" )->SetDescription( "Mask file name" )->SetRequired( true )->SetMinMax( 1, 1 );
	p.AddArgument( list.outputFileName, "output" )->AddAlias( "o" )->SetDescription( "Output file name" )->SetRequired( true )->SetMinMax(
			1, 1 );
	p.AddArgument( list.threshold, "threshold" )->AddAlias( "t" )->SetDescription( "Fixed degree threshold (default: 50)" )->SetRequired(
			false )->SetMinMax( 1, 1 );

	p.AddArgument( list.noThreshold, "no-threshold" )->AddAlias( "nt" )->SetDescription( "No threshold" )->SetRequired(
				false )->SetMinMax( 1, 1 );

	p.AddArgument( list.clusters, "clusters" )->AddAlias( "c" )->SetDescription( "Number of clusters (default: 5)" )->SetRequired( false )->SetMinMax(
			1, 1 );
	p.AddArgument( list.repeat, "repeat" )->AddAlias( "r" )->SetDescription( "Repeat clustering (default: 10)" )->SetRequired( false )->SetMinMax(
			1, 1 );
	p.AddArgument( list.verbose, "verbose" )->AddAlias( "v" )->SetDescription( "Print status to screen (default: false)" )->SetRequired(
			false )->SetMinMax( 1, 1 );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnClustering annClustering;
	annClustering.Run( list );

	return EXIT_SUCCESS;
}

