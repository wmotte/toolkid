#include "tkdCmdParser.h"

#include "annCommon.h"

/**
 * Calculate modularity.
 */
class AnnModularity
{
public:

	typedef float PixelType;
	typedef std::vector< PixelType > VectorType;
	typedef std::vector< VectorType > MatrixType;

	typedef itk::Image< PixelType, 2 > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
	typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIteratorType;

	typedef boost::tokenizer< boost::char_separator< char > > TokType;

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName, const std::string& clusterFileName )
	{
		// if input is 2D image...
		if ( ann::Ann< PixelType >::GetImageDimensions( inputFileName ) == 2 )
		{
			ann::Ann< PixelType >::MatrixType matrix;
			ann::Ann< PixelType >::LoadMatrix( matrix, inputFileName );
			ann::Ann< PixelType >::VectorType clusters;
			ann::Ann< PixelType >::LoadVector( clusters, clusterFileName );
			std::cout << ann::Ann< PixelType >::Modularity( clusters, matrix ) << std::endl;
		}
		else
		{
			std::cerr << "*** ERROR ***: Dimensions of input image is not equal to 2 (matrix)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		exit( EXIT_SUCCESS );
	}
};

/**
 * Modularity.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Calculate clustering modularity" );

	std::string inputFileName;
	std::string clusterFileName;

	p.AddArgument( inputFileName, "input" )->AddAlias( "i" )->SetDescription( "Input file name (2D matrix)" )->SetRequired( true )->SetMinMax(
			1, 1 );
	p.AddArgument( clusterFileName, "cluster" )->AddAlias( "c" )->SetDescription( "Cluster file name (txt)" )->SetRequired( true )->SetMinMax(
			1, 1 );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnModularity annModularity;
	annModularity.Run( inputFileName, clusterFileName );

	return EXIT_SUCCESS;
}

