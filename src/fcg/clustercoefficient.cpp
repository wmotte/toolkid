#include "itkImage.h"
#include "itkImageFileReader.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "tkdCmdParser.h"
#include "flens.h"

class ClusterCoefficient {
public:
	typedef float PixelType;
	typedef itk::Image< PixelType, 2 > ImageType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef float PrecisionType;
	typedef flens::GeMatrix< flens::FullStorage< PrecisionType, flens::RowMajor > > MatrixType;
	typedef vnl_vector< PrecisionType > VectorType;

	void RunWeighted( const std::string& filename, PixelType threshold ) {
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();
		reader = 0;

		PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
		ImageType::RegionType region = image->GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		int rows = size[0];
		int cols = size[1];

		DataMatrixType wIn( rows, cols, buffer );

		MatrixType A( rows, cols );
		MatrixType S( rows, cols );
		MatrixType B( rows, cols );
		MatrixType W( rows, cols );

		for ( int i = 0; i < rows; ++i ) {
			for ( int j = 0; j < cols; ++j ) {
				S( i + 1, j + 1 ) = 0;

				if ( wIn( i, j ) > threshold ) {
					A( i + 1, j + 1 ) = 1;
					B( i + 1, j + 1 ) = vcl_pow( wIn( i, j ), 1. / 3. );
					W( i + 1, j + 1 ) = wIn( i, j );
				} else {
					A( i + 1, j + 1 ) = 0;
					B( i + 1, j + 1 ) = 0;
					W( i + 1, j + 1 ) = 0;
				}
			}
		}

		for ( int i = 1; i <= rows; ++i ) {
			for ( int j = 1; j <= cols; ++j ) {
				S( i, j ) = B( i, j ) + B( j, i );
			}
		}

		MatrixType p1( rows, cols );
		MatrixType p( rows, cols );
		MatrixType a( rows, cols );
		VectorType C( rows );
		VectorType K( rows );

		flens::copy( S * S, p1 );
		flens::copy( p1 * S, p );
		flens::copy( A * A, a );

		K.fill( 0 );
		for ( int i = 0; i < rows; ++i ) {
			for ( int j = 0; j < cols; ++j ) {
				K( i ) += ( A( i + 1, j + 1 ) + A( j + 1, i + 1 ) );
			}

			PrecisionType cyc = p( i + 1, i + 1 ) / 2.;
			if ( cyc == 0 ) {
				K( i) = vcl_numeric_limits< PrecisionType>::max();
			}

			C( i ) = cyc / ( K( i ) * ( K( i ) - 1 ) - 2 * a( i + 1, i + 1) );
		}

		PrecisionType sum = 0;
		for( int i = 0; i < rows; ++i )
		{
			sum += C( i );
		}

		sum /= static_cast< PrecisionType>( rows );
		std::cout << sum << std::endl;
	}

	void Run( const std::string& filename, PixelType threshold )
	{
		typedef itk::ImageFileReader< ImageType> ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();
		reader = 0;

		PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
		ImageType::RegionType region = image->GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		int rows = size[ 0 ];
		int cols = size[ 1 ];

		DataMatrixType data( rows, cols, buffer );
		for( int i = 0; i < ( rows * cols ); ++i )
		{
			buffer[ i ] *= ( buffer[ i ] < 0 ? -1. : 1. );
		}

		PixelType sum = 0;
		//
						//    std::vector< PixelType > t;
						//    for( int i = 0; i < rows; ++i )
						//      {
						//      for( int j = i + 1; j < cols; ++j )
						//        {
						//        if ( data( i, j ) > 0 )
						//          {
						//          t.push_back( data( i, j ) );
						//          }
						//        }
						//      }
						//
						//    std::sort( t.begin(), t.end() );
						//    int index = static_cast< int >( 0.8 * static_cast< float >( t.size() ) );
						//    std::cout << t[ index ] << std::endl;
						//    return;

						for( int i = 0; i < rows; ++i )
						{
							// 1) determine k = number of edges to neighbors
							int k = 0;
							std::vector< int> neighbors;
							for( int j = 0; j < cols; ++j )
							{
								if ( i == j )
								{
									continue;
								}

								if ( data( i, j ) < threshold )
								{
									continue;
								}

								++k;
								neighbors.push_back( j );
							}

							// 2) determine connections among neighbors
							int connections = 0;
							std::vector< int>::const_iterator end = neighbors.end();
							for( std::vector< int>::const_iterator j = neighbors.begin(); j != end; ++j )
							{
								for( std::vector< int>::const_iterator m = j; m != end; ++m )
								{
									if ( *j == *m || i == *j || *m == i )
									{
										continue;
									}

									if ( data( *j, *m ) < threshold )
									{
										continue;
									}

									++connections;
								}
							}

							// k(k-1)*0.5 possible connections
							if ( k> 1 )
							{
								PixelType coefficient = static_cast< PixelType>( 2 * connections ) / static_cast< PixelType>( k * ( k - 1 ) );
								sum += coefficient;
								//        std::cout << i << "\t" << coefficient << std::endl;
							}
						}

						PixelType average = sum / static_cast< PixelType>( rows );

						std::cout << average << std::endl;
					}
				};

int main( int argc, char ** argv ) {
	tkd::CmdParser p( "clustercoefficient", "Calculate voxel-wise cluster coefficients" );

	std::string inputFileName;
	float threshold;
	bool weightedGraph = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Input 2D image: adjacency matrix" ) ->SetRequired(
			true );

	p.AddArgument( threshold, "threshold" ) ->AddAlias( "t" ) ->SetInput( "float" ) ->SetDescription( "Threshold" ) ->SetRequired( true );

	p.AddArgument( weightedGraph, "weighted" ) ->AddAlias( "w" ) ->SetDescription( "Calculated weighted cluster coefficients" );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return -1;
	}

	ClusterCoefficient cc;

	if ( weightedGraph ) {
		cc.RunWeighted( inputFileName, threshold );
	} else {
		cc.Run( inputFileName, threshold );
	}

	return 0;
}
