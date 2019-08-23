#include "tkdCmdParser.h"

#include "nrFourier.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_fft_1d.h"
#include "vnl/algo/vnl_cholesky.h"

#include "itkNumericTraits.h"

#include <stdlib.h>
#include <math.h>
#include <fstream>

#include <boost/math/constants/constants.hpp>

#include "Array.h"
#include "fftw++.h"

using namespace Array;
using namespace fftwpp;

/**
 * Option list.
 */
struct parameters
{
	std::string inputFileName;
	std::string separation;
	bool verbose;
};

/**
 * General Linear Model (GLM) with autocorrelation prewhitening.
 *
 * Refs:
 * 1. Book: 'Functional MRI an introduction to methods',
 * 		P. Jezzard, P.M. Matthews, S.M. Smith, Oxford (2005).
 *
 * 2. Article: 'Temporal autocorrelation in Univariate Modeling of FMRI Data',
 * 		M.W. Woolrich, B.D. Ripley, M. Brady and S.M. Smith, Neuroimage 14, 1370-1386 (2001).
 */
namespace general_linear_model
{
	typedef double PixelType;
	typedef vnl_matrix< PixelType > MatrixType;
	typedef vnl_vector< PixelType > VectorType;

	class GLM
	{
	public:

		/**
		 * Run glm.
		 */
		void Run( parameters& list )
		{
			// data Y
			VectorType Y( 8, 0 );
			Y( 0 ) = 0.1;
			Y( 1 ) = 3.5;
			Y( 2 ) = 6.7;
			Y( 3 ) = 6.4;
			Y( 4 ) = 3.5;
			Y( 5 ) = 3.7;
			Y( 6 ) = 3.1;
			Y( 7 ) = 2.1;

			// design matrix X
			MatrixType X( 8, 3 );
			X( 0, 0 ) = 1;
			X( 1, 0 ) = 1;
			X( 2, 0 ) = 1;
			X( 3, 0 ) = 1;
			X( 4, 0 ) = 1;
			X( 5, 0 ) = 1;
			X( 6, 0 ) = 1;
			X( 7, 0 ) = 1;

			X( 0, 1 ) = 0;
			X( 1, 1 ) = 0;
			X( 2, 1 ) = 1;
			X( 3, 1 ) = 1;
			X( 4, 1 ) = 0;
			X( 5, 1 ) = 0;
			X( 6, 1 ) = 0;
			X( 7, 1 ) = 0;

			X( 0, 2 ) = 0;
			X( 1, 2 ) = 0;
			X( 2, 2 ) = 0;
			X( 3, 2 ) = 0;
			X( 4, 2 ) = 1;
			X( 5, 2 ) = 1;
			X( 6, 2 ) = 0;
			X( 7, 2 ) = 0;

			// Toeplitz matrix S
			MatrixType S( Y.size(), Y.size() );
			S.set_identity();

			// Adjust for correlation using Cochrane–Orcutt estimation
			VectorType B = OLS( Y, X, S );

			std::cout << "Betas: " << B << std::endl;

			// Residuals
			VectorType r = Residuals( Y, X, B );

			std::cout << "Residuals: " << r << std::endl;

			// Raw autocorrelation (all lags)
			VectorType ac = AutoCorrelation( r );

			std::cout << "Raw autocorrelation: " << ac << std::endl;

			// non-parametric autocorrelation smoothing (M = 2sqrt(N))
			unsigned int M = 2. * std::sqrt( Y.size() );

			VectorType sac = TaperedCosineWindow( ac, M );

			std::cout << "Smoothed autocorrelation (" << M << ")" << ": " << sac << std::endl;

			// construct V
			MatrixType V = GetV( sac, true );

			vnl_svd< PixelType > svd( V );
			//svd.solve();
			MatrixType W = svd.W();
			MatrixType Vprime = W * W.transpose();



			for( unsigned int r = 0; r < V.rows(); r++ )
			{
				for( unsigned int c = 0; c < V.cols(); c++ )
				{
					std::cout << V( r, c ) << " ";
				}
				std::cout << std::endl;
			}

			std::cout << "UU'" << std::endl;

			for( unsigned int r = 0; r < Vprime.rows(); r++ )
			{
				for( unsigned int c = 0; c < Vprime.cols(); c++ )
				{
					std::cout << Vprime( r, c ) << " ";
				}
				std::cout << std::endl;
			}


			exit( 0 );

			// refit
			//B = OLS( Y, X, S );
		}

	protected:

		/**
		 * Return inverted Cholesky V=KK'
		 */
		MatrixType GetInvertedK( const MatrixType& V )
		{
			vnl_svd< PixelType > svd( V );
			MatrixType pinvV = svd.inverse();


		}

		/**
		 * Return V.
		 */
		MatrixType GetV( const VectorType& v, bool circular )
		{
			MatrixType result( v.size(), v.size(), 0.0 );

			for( unsigned int r = 0; r < result.rows(); r++ )
				for( unsigned int c = r; c < result.cols(); c++ )
					result( r, c ) = v( r );

			result += result.transpose();

			result.fill_diagonal( 1 );

			return result;
		}


		/**
		 * Tukey windowing.
		 *
		 * http://en.wikipedia.org/wiki/Window_function#Tukey_window
		 */
		VectorType TaperedCosineWindow( const VectorType& ar, unsigned int M )
		{
			VectorType rho( ar.size(), 0 );

			for ( unsigned int i = 0; i < ar.size(); i++ )
			{
				if ( i < M )
				{
					PixelType PI = boost::math::constants::pi< PixelType >();
					rho( i ) = 0.5 * ( 1 + std::cos( ( PI * i ) / M ) ) * ar( i );
				} else
				{
					rho( i ) = 0;
				}
			}

			return rho;
		}

		/**
		 * Cross correlation for all lags
		 * @param a Vector 1 (of length n)
		 * @param b Vector 2 (of length n)
		 * @param output Output vector of length n
		 *
		 * Normalize the sequence so the autocorrelations at zero lag are identically 1.0.
		 */
		VectorType AutoCorrelation( const VectorType& a )
		{
			const int n = a.size();

			// added zeropadding!
			int zeropad = (int) pow( 2, ceil( log( a.size() ) / log( 2 ) ) );
			VectorType output( zeropad, 0 );

			for ( int i = 0; i < n; ++i )
				output[i] = a[i];

			array1< PixelType > finput( a.size(), sizeof( Complex ) );
			array1< Complex > foutput( a.size(), sizeof( Complex ) );
			array1< Complex > ftemp( a.size(), sizeof( Complex ) );

			for( unsigned int i = 0; i < a.size(); i++ )
				finput[ i ] = a( i );

			rcfft1d Forward( a.size(), finput, foutput );

			Forward.fft( finput, foutput ); // fill foutput ...
			Forward.fft( finput, ftemp ); // fill ftemp ...

			const PixelType no2 = static_cast< PixelType > ( n >> 1 );

			for ( int i = 2; i < n; i += 2 )
			{
				Complex ftmp = foutput[i];
				foutput[i] = ( foutput[i] * ftemp[i] + foutput[i + 1] * ftemp[i + 1] ) / no2;
				foutput[i + 1] = ( foutput[i + 1] * ftemp[i] - ftmp * ftemp[i + 1] ) / no2;
			}

			foutput[0] = foutput[0] * ftemp[0] / no2;
			foutput[1] = foutput[1] * ftemp[1] / no2;

			crfft1d Backward( a.size(), foutput, finput );
			Backward.fftNormalized( foutput, finput );

			VectorType finalOutput( n );

			// recrop
			for ( int i = 0; i < n; ++i )
				finalOutput( i ) = std::abs( foutput( i ) / foutput( 0 ) ); // normalize

			return finalOutput;
		}

		/**
		 * Return standard deviation.
		 */
		template< class T >
		T GetSD( const vnl_vector< T >& distances )
		{
			int size = distances.size();

			if ( size < 2 )
			{
				return 0;
			}

			T mean = distances.mean();

			T sd = 0;
			for ( int i = 0; i < size; i++ )
			{
				sd += pow( distances[i] - mean, 2 );
			}

			return vcl_sqrt( sd / size );
		}

		/**
		 * Return residuals of fit.
		 */
		VectorType Residuals( const VectorType& Y, const MatrixType& X, const VectorType& B )
		{
			return Y - X * B;

			// TODO
			//VectorType r( B.size(), 0 );
			//return r;
		}

		/**
		 * Ordinary Least Squares ( OLS ) estimation, by finding the
		 * maximum likelihood estimates of a linear regression model.
		 *
		 * Ref: http://en.wikipedia.org/wiki/Ordinary_least_squares
		 */
		VectorType OLS( const VectorType& Y, const MatrixType& X, const MatrixType& S )
		{
			if ( Y.size() != X.rows() )
			{
				std::cerr << "*** ERROR ***: The data vector Y and rows of the "
					"design matrix X must have the same length." << std::endl;
				exit( EXIT_FAILURE );
			}
			if ( X.rows() != S.rows() )
			{
				std::cerr << "*** ERROR ***: The Toeplitz matrix S and "
					"design matrix X must have the same number of rows." << std::endl;
				exit( EXIT_FAILURE );
			}

			// SY = SXB + mu, where mu ~ N( 0, sigma^2 SVS') => Bhat = pinv(SX)SY,
			// where pinv(SX) is the Moore-Penrose pseudoinverse: ((SX)'SX)^-1(SX)'
			vnl_svd< PixelType > svd( ( S * X ).transpose() * S * X );
			MatrixType pinvSX = svd.inverse() * ( S * X ).transpose();

			return pinvSX * S * Y;
		}
	};
} // end namespace general_linear_model


/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Generalized linear model" );

	parameters list;
	list.verbose = false;

	p.AddArgument( list.inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription(
			"Input file: list of covariate(s) and response function (column format)" ) ->SetRequired( true );

	p.AddArgument( list.separation, "separation" ) ->AddAlias( "s" ) ->SetInput( "string" ) ->SetDescription(
			"Data separator string (default: \\t" );

	p.AddArgument( list.verbose, "verbose" ) ->AddAlias( "v" ) ->SetInput( "bool" ) ->SetDescription( "Verbose (default: false" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	general_linear_model::GLM glm;

	glm.Run( list );

	return EXIT_SUCCESS;
}

