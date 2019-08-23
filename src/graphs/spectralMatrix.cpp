#include "spectralMatrix.h"

namespace fcg
{
	/**
	 * Run.
	 */
	void SpectralMatrix::Run( const std::string& inputFileName, const std::string& maskFileName,
			const std::string& outputFileName,
			int skipBegin, int skipEnd, double sigma, bool negativeEuclidean )
	{
		std::cout << "Read "<< inputFileName << std::endl;
		ImageType::Pointer image = LoadImage( inputFileName );
		std::cout << "Read "<< maskFileName << std::endl;
		MaskType::Pointer mask = LoadMask( maskFileName );

		ImageType::RegionType region = image->GetLargestPossibleRegion( );
		ImageType::SizeType size = region.GetSize( );

		int voxels = size[ 0 ] * size[ 1 ]* size[ 2 ];
		int timepoints = size[ 3 ];

		PixelType* data = image->GetPixelContainer()->GetBufferPointer( );
		MaskPixelType* dataMask = const_cast< MaskType::ImageType* >( mask->GetImage() )->GetPixelContainer()->GetBufferPointer( );

		DataType source( timepoints, voxels, data );
		MatrixType matrix = source.extract( source.rows( ) - skipBegin - skipEnd, source.cols( ), skipBegin, 0 );

		std::vector< int > indices;
		std::vector< int > inverseIndices;
		int rowCount = 0;
		for ( int i = 0; i < voxels; ++i )
		{
			if ( dataMask[ i ] )
			{
				indices.push_back( i );
				inverseIndices.push_back( rowCount++ );
			} else
			{
				inverseIndices.push_back( -1 );
			}
		}

		int validVoxels = indices.size();
		std::cout << validVoxels << " voxels with "<< matrix.rows( ) << " timepoints"<< std::endl;

		OutputImageType::Pointer outputImage = OutputImageType::New( );
		OutputImageType::RegionType outputRegion;
		OutputImageType::SizeType outputSize;
		outputSize[ 0 ] = validVoxels;
		outputSize[ 1 ] = validVoxels;
		outputRegion.SetSize( outputSize );
		outputImage->SetRegions( outputRegion );
		outputImage->Allocate( );
		outputImage->FillBuffer( itk::NumericTraits<OutputPixelType>::Zero );
		OutputPixelType* outputData = outputImage->GetPixelContainer()->GetBufferPointer( );
		OutputDataType output( validVoxels, validVoxels, outputData );

		std::vector< VectorType > vectors;

		for ( int i = 0; i < validVoxels; ++i )
		{
			VectorType a = matrix.get_column( indices[ i ] );
			vectors.push_back( a );
		}

		image = 0;

		std::cout << "Calculate functional distance" << std::endl;

		int numberOfCalculations = ( validVoxels ) * validVoxels * 0.5;
		int counter = 0;
		int lastPercentage = 0;
		std::vector< OutputPixelType > medianVector;

		// first, calculate functional distance ...
		for( int i = 0; i < validVoxels; ++i )
		{
			const VectorType& a = vectors[ i ];

			for( int j = i; j < validVoxels; ++j )
			{
				const VectorType& b = vectors[ j ];

				OutputPixelType r = 0;

				if( !negativeEuclidean )
				{
					// get raw functional distance ...
					r = nr::correlation::FunctionalDistance< OutputPixelType >( a, b );
				}
				else
				{
					// negative euclidean distance ...
					r = - nr::correlation::FunctionalDistance< OutputPixelType >( a, b );
				}

				medianVector.push_back( r );

				output( i, j ) = r;
				output( j, i ) = r;

				int percentage = ( static_cast< double >( ++counter ) / static_cast< double >( numberOfCalculations ) ) * 100.;
				if ( percentage > lastPercentage )
				{
					std::cout << "\r" << percentage << "%";
					std::cout.flush();
					lastPercentage = percentage;
				}
			}
		}
		std::cout << std::endl;

		// set sigma to median if sigma is nog manual defined...
		if( sigma > std::numeric_limits< double >::min() )
		{
			sigma = graph::Graph< OutputPixelType >::GetMedian( medianVector );
			std::cout << "Median: " << sigma << std::endl;
		}

		// third, normalize using scaling factor ...
		std::cout << "Normalizing" << std::endl;
		lastPercentage = 0;

		// normalize if sigma is manually defined...
		if( sigma > std::numeric_limits< double >::min() )
		{
			for( int i = 0; i < validVoxels; ++i )
			{
				for( int j = i; j < validVoxels; ++j )
				{
					OutputPixelType norm = std::exp( - std::pow( output( i, j ) / sigma, static_cast< OutputPixelType >( 2 ) ) );

					output( i, j ) = norm;
					output( j, i ) = norm;
				}
			}
		}

		std::cout << std::endl;

		typedef itk::ImageFileWriter< OutputImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( outputFileName.c_str() );
		writer->SetInput( outputImage );

		std::cout << "Write " << outputFileName << std::endl;
		writer->Update();
	}

	/**
	 * Load image.
	 */
	SpectralMatrix::ImageType::Pointer SpectralMatrix::LoadImage( const std::string& filename )
	{
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		return reader->GetOutput();
	}

	/**
	 * Load mask.
	 */
	SpectralMatrix::MaskType::Pointer SpectralMatrix::LoadMask( const std::string& filename )
	{
		typedef itk::ImageFileReader< MaskImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		MaskType::Pointer mask = MaskType::New();
		mask->SetImage( reader->GetOutput() );
		return mask;
	}
} // end namespace fcg

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "spectral_matrix", "Calculate similarity matrix for spectral clustering");

	std::string inputFileName;
	std::string maskFileName;
	std::string outputFileName;
	double sigma = std::numeric_limits< double >::min();
	bool negativeEuclidean = false;

	p.AddArgument( inputFileName, "input" )
	->AddAlias( "i" )
	->SetInput( "filename" )
	->SetDescription( "Input 4D time series image" )
	->SetRequired( true );

	p.AddArgument( maskFileName, "mask" )
	->AddAlias( "m" )
	->SetInput( "filename" )
	->SetDescription( "Input 3D mask image" )
	->SetRequired( true );

	p.AddArgument( outputFileName, "output" )
	->AddAlias( "o" )
	->SetInput( "filename" )
	->SetDescription( "Output 2D image: similarity matrix" )
	->SetRequired( true );

	p.AddArgument( sigma, "sigma" )
	->SetInput( "float" )
	->SetDescription( "Scaling factor in functional distance measure (default: median, if no non-zero value is specified)" )
	->SetRequired( false );

	p.AddArgument( negativeEuclidean, "negative-euclidean" )
	->SetInput( "bool" )
	->SetDescription( "Fill matrix with negative euclidean distance (default: false)" )
	->SetRequired( false );

	std::vector< int > skip;
	skip.push_back( 50 );
	skip.push_back( 50 );

	p.AddArgument( skip, "skip" )
	->AddAlias( "s" )
	->SetMinMax( 2, 2 )
	->SetDescription( "Skip samples at begin/end (default: 50 50)" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	fcg::SpectralMatrix spectralMatrix;
	spectralMatrix.Run( inputFileName, maskFileName, outputFileName, skip[ 0 ], skip[ 1 ], sigma, negativeEuclidean );

	return 0;
}
