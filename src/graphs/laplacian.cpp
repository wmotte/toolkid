#include "tkdCmdParser.h"

#include "graphCommon.h"

/**
 * Calculate (normalized) laplacian and eigenvector/values.
 */
class Laplacian
{
public:

	typedef double PixelType;

	/**
	 * Run app.
	 */
	void Run( const std::string& filename, const std::string& outputFileName,
			bool normalize,
			const std::string& outputFileNameEigen,
			unsigned int numberOfEigenVectors,
			bool outputEigenFiles,
			const std::string& outputFileNameKmeans,
			unsigned int numberOfClusters,
			bool outputKmeans,
			bool sparse )
	{
		// get dimensions...
		unsigned int dims = graph::Graph< PixelType >::GetImageDimensions( filename );

		if ( dims == 2 )
		{
			graph::Graph< PixelType >::ImagePointerType image;
			graph::Graph< PixelType >::ReadMatrix( image, filename );

			PixelType* buffer = image->GetPixelContainer()->GetBufferPointer( );
			graph::Graph< PixelType >::ImageType::RegionType region = image->GetLargestPossibleRegion( );
			graph::Graph< PixelType >::ImageType::SizeType size = region.GetSize();
			int rows = size[ 0 ];
			int cols = size[ 1 ];

			graph::Graph< PixelType >::DataMatrixType G( rows, cols, buffer );
			graph::Graph< PixelType >::DiagMatrixType D;
			graph::Graph< PixelType >::MatrixType L;

			// [ 1 ]
			graph::Graph< PixelType >::CalculateLaplacian( L, G, D, normalize );

			// [ 2 ]
			graph::Graph< PixelType >::MatrixType eigenVectors;
			graph::Graph< PixelType >::VectorType eigenValues;

			if( sparse )
			{
				// convert dense matrix into sparse matrix...
				vnl_sparse_matrix< PixelType > S( L.rows(), L.cols() );
				for( unsigned int i = 0; i < L.rows(); i++ )
					for( unsigned int j = i; j < L.cols(); j++ )
						if( L( i, j ) != 0 )
							S( i, j ) = L( i, j );

				graph::Graph< PixelType >::GeneralizedEigenCalculation( S, eigenVectors, eigenValues, numberOfEigenVectors );
			}
			else
				graph::Graph< PixelType >::GeneralizedEigenCalculation( L, D, eigenVectors, eigenValues );


			// Output eigenfiles...
			graph::Graph< PixelType >::WriteEigenFiles( outputFileNameEigen, numberOfEigenVectors, eigenVectors, eigenValues );

			if( !outputFileName.empty() )
			{
				graph::Graph< PixelType >::WriteMatrixToFile( L, outputFileName );
			}

			// [ 3 ] Optional: run KMEANS++
			if( outputKmeans )
			{
				graph::Graph< PixelType >::VectorType clusters;
				graph::Graph< PixelType >::KMeansPP( eigenVectors, numberOfClusters, clusters );

				graph::Graph< PixelType >::WriteVectorToFile(
						outputFileNameKmeans + "_clusters.txt", clusters );
			}
		}
		else
		{
			std::cerr << "Number of input (matrix) dimensions should be 2!"<< std::endl;
			exit( EXIT_FAILURE );
		}
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Compute Laplacian of graph's adjacency matrix");

	std::string inputFileName;
	std::string outputFileName;

	std::string outputFileNameEigen = "eigen";
	int numberOfEigenVectors = 15;
	bool outputEigenFiles = false;

	bool normalize = false;

	std::string outputFileNameKmeans = "kmeans";
	int numberOfClusters = 15;
	bool outputKmeans = false;
	bool sparse = false;

	p.AddArgument( inputFileName, "input" )
	->AddAlias( "i" )
	->SetInput( "filename" )
	->SetDescription( "Input image: 2D adjacency matrix" )
	->SetRequired( true );

	p.AddArgument( outputFileName, "output" )
	->AddAlias( "o" )
	->SetInput( "filename" )
	->SetDescription( "Output image: 2D laplacian matrix" )
	->SetRequired( false );

	p.AddArgument( outputFileNameEigen, "output-eigen-root" )
	->AddAlias( "oer" )
	->SetInput( "root filename" )
	->SetDescription( "Root eigen files (default: 'eigen')" )
	->SetRequired( false );

	p.AddArgument( numberOfEigenVectors, "output-eigenvector-number" )
	->AddAlias( "oen" )
	->SetInput( "int" )
	->SetDescription( "Number of eigenvectors to output (default: 15; first (zero) eigenvector is ignored)" )
	->SetRequired( false );

	p.AddArgument( outputEigenFiles, "output-eigen-files" )
	->AddAlias( "oef" )
	->SetInput( "bool" )
	->SetDescription( "Output eigenvectors and eigenvalues (default: false)" )
	->SetRequired( false );

	p.AddArgument( outputFileNameKmeans, "output-kmeans-root" )
	->AddAlias( "okr" )
	->SetInput( "root filename" )
	->SetDescription( "Root kmeans files (default: 'kmeans')" )
	->SetRequired( false );

	p.AddArgument( numberOfClusters, "output-kmeans-clusters" )
	->AddAlias( "okc" )
	->SetInput( "int" )
	->SetDescription( "Number of clusters to output (default: 15)" )
	->SetRequired( false );

	p.AddArgument( outputKmeans, "output-kmeans-file" )
	->AddAlias( "okf" )
	->SetInput( "bool" )
	->SetDescription( "Output kmeans++ segmentation (default: false)" )
	->SetRequired( false );

	p.AddArgument( normalize, "normalize" )
	->AddAlias( "n" )
	->SetInput( "bool" )
	->SetDescription( "Normalize laplacian matrix (default: false)" )
	->SetRequired( false );

	p.AddArgument( sparse, "sparse" )
	->AddAlias( "s" )
	->SetInput( "bool" )
	->SetDescription( "Use Lanczos algorithm to estimate eigen-vectors (default: false)" )
	->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Laplacian laplacian;
	laplacian.Run( inputFileName, outputFileName, normalize,
			outputFileNameEigen, (unsigned) numberOfEigenVectors, outputEigenFiles,
			outputFileNameKmeans, (unsigned) numberOfClusters, outputKmeans, sparse );

	return EXIT_SUCCESS;
}
