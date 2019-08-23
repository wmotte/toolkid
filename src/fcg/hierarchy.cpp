#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/floyd_warshall_shortest.hpp"
#include "boost/graph/adjacency_matrix.hpp"

#include <vnl/algo/vnl_lsqr.h>
#include <vnl/vnl_sparse_matrix_linear_system.h>
#include <vnl/vnl_least_squares_function.h>
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"

#include "tkdCmdParser.h"

#include "flens.h"

#include <iostream>

/**
 * Beta parameter of whole-weighted matrix is estimated, according to:
 *
 * Ravasz E, Barabasi AL (2003) Hierarchical organization in complex networks.
 * Phys Rev E Stat Nonlin Soft Matter Phys 67: 026112
 *
 * Beta measures the extent of the power-law relationship between the clustering coefficient (C) and
 * the degree (k): C ~ k^-beta.
 *
 * We calculate the weighted clustering coefficient and weighted degree ( ~ strength ) for each
 * voxel and fit a linear regression line to the plot of log(C) versus log(k).
 *
 * Date: 29-09-2009
 */
class Hierarchy {

public:

	typedef float PixelType;
	typedef double PrecisionType;

	/**
	 * Run.
	 */
	void run( const std::string& inputFileName, PixelType threshold, const std::string& csvFileName ) {

		// get dimensions...
		unsigned int dims = getDimensions( inputFileName );

		if ( dims == 2 ) {

			// [ 1 ] get weighted clustering-coefficients
			std::vector< PrecisionType > ccs = ClusteringCoefficients( inputFileName, threshold );

			// [ 2 ] get degrees
			std::vector< PrecisionType > degrees = Degrees( inputFileName, threshold );

			// [ 3 ] least squares fit log(C) versus log(k)
			fit( ccs, degrees, csvFileName );

		} else {
			std::cerr << "Number of input (matrix) dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	typedef itk::Image< PixelType, 2 > ImageType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef flens::GeMatrix< flens::FullStorage< PrecisionType, flens::RowMajor > > MatrixType;
	typedef vnl_vector< PrecisionType > VectorType;

	/**
	 * Check dimensions of inputfile. In case of error, EXIT_FAILURE is returned.
	 */
	unsigned int getDimensions( const std::string& inputFileName ) {
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io ) {
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else {
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Return degrees.
	 */
	std::vector< PrecisionType > Degrees( const std::string& filename, PixelType threshold ) {

		typedef itk::Image< PixelType, 2 > ImageType;
		typedef vnl_matrix_ref< PixelType > DataMatrixType;

		typedef flens::GeMatrix< flens::FullStorage< PrecisionType, flens::RowMajor > > MatrixType;
		typedef vnl_vector< PrecisionType > VectorType;

		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( filename );
		reader -> Update();

		ImageType::Pointer image = reader -> GetOutput();
		reader = 0;

		PixelType* buffer = image -> GetPixelContainer() -> GetBufferPointer();
		ImageType::RegionType region = image -> GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		int rows = size[0];
		int cols = size[1];

		DataMatrixType data( rows, cols, buffer );
		for ( int i = 0; i < ( rows * cols ); ++i ) {
			buffer[i] *= ( buffer[i] < 0 ? -1. : 1. );
		}

		// fill rows with sum of column-values...
		std::vector< PrecisionType > degrees( rows, 0 );

		for ( int i = 0; i < rows; i++ ) {
			degrees[i] = 0;
			for ( int j = 0; j < cols; j++ ) {
				degrees[i] += data( i, j );
			}
		}
		return degrees;
	}

	/**
	 * Least squares fit log(C) versus log(k).
	 *
	 * fit linear regression line to the plot of log(C) versus log(k).
	 * and output results.
	 */
	void fit( const std::vector< PrecisionType >& clusteringCoefficients, const std::vector< PrecisionType > degrees, const std::string& csvFileName ) {

		// sanity check...
		if ( degrees.size() != clusteringCoefficients.size() ) {
			std::cout << "Size of degrees does not match size of clustering coefficients!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// only fit non-thresholded voxels...
		unsigned int N = 0;
		for( unsigned int i = 0; i < clusteringCoefficients.size(); i++ ) {
			if( clusteringCoefficients[i] > 0 ) {
				N++;
			}
		}

		// init...
		vnl_sparse_matrix< PrecisionType > A( N, 1 );
		vnl_vector< PrecisionType > b( N );

		// prepare...
		unsigned int count = 0;
		for ( unsigned int i = 0; i < clusteringCoefficients.size(); i++ ) {

			double x = degrees[i];
			double y = clusteringCoefficients[i];

			if ( x > 0 ) {
				A( count, 0 ) = log( x );
				b[count]      = log( y );
				count++;
			}
		}

		vnl_sparse_matrix_linear_system< PrecisionType > ls( A, b );
		vnl_vector< PrecisionType > x( 1 );
		x[0] = 0.0;
		vnl_lsqr lsqr( ls );

		// fit...
		lsqr.minimize( x );

		double residuals = lsqr.get_resid_norm_estimate();
    double beta = -1 * x[0];
		std::cout << "beta: " << beta << std::endl;
		std::cout << "residual: " << residuals << std::endl;

		// write raw data to txt-file
		if ( csvFileName.size() != 0 ) {

			std::ofstream filestream;
			filestream.open( csvFileName.c_str() );
			if ( filestream.fail() ) {
				std::cerr << "Not able to write to: " << csvFileName << std::endl;
				exit( EXIT_FAILURE );
			}
			else {

				filestream << "Lsqr beta: " << beta << std::endl;
				filestream << "Lsqr res.: " << residuals << std::endl << std::endl;
				filestream << "Voxel, Weighted clustering coefficient,Weighted degree" << std::endl;

				for ( unsigned int i = 0; i < clusteringCoefficients.size(); i++ ) {

					double C = clusteringCoefficients[i];
					double k = degrees[i];

					if ( C > 0 ) {
						filestream << i << "," << C << "," << k << std::endl;
					}
				}
			}

			filestream.close();
		}
	}

	/**
	 * Return clustering coefficients.
	 */
	std::vector< PrecisionType > ClusteringCoefficients( const std::string& filename, PixelType threshold ) {

		std::vector< PrecisionType > clusteringCoefficients;

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

		for( int i = 0; i < rows; ++i )
		{
			clusteringCoefficients.push_back( C( i ) );
		}

		return clusteringCoefficients;
	}
};

/**
 * Main hierarchy.
 */
int main( int argc, char ** argv ) {
	tkd::CmdParser p( "hierarchy", "Calculate hierarchy parameter beta (weighted graphs only)" );

	std::string inputFileName;
	std::string csvFileName;
	float threshold = 0;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Input 2D image: adjacency matrix" ) ->SetRequired(
			true );

	p.AddArgument( csvFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output raw data" ) ->SetRequired(
			false );

	p.AddArgument( threshold, "threshold" ) ->AddAlias( "t" ) ->SetInput( "float" ) ->SetDescription( "Threshold [ default: 0 ]" ) ->SetRequired( false );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Hierarchy hierarchy;

	hierarchy.run( inputFileName, threshold, csvFileName );

	return EXIT_SUCCESS;
}

