#include "tkdCmdParser.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matrix_fixed.h"
#include <vnl/algo/vnl_cholesky.h>
#include <vnl/vnl_transpose.h>

#include <cmath>
#include <ctime>

#include <list>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkDiffusionTensor3D.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/math/distributions/gamma.hpp>

/**
 * dki namespace.
 */
namespace dki
{

	struct parameters
	{
		std::string inputFileName;
		std::string maskFileName;
		std::string bvecsFileName;
		std::string bvalsFileName;
		std::string outputFileName;
		std::string sep;
		int iterations;
		bool regularize;
	};

	typedef double PixelType;
	typedef unsigned char MaskPixelType;

	typedef itk::Image< PixelType, 4 > ImageType;
	typedef itk::Image< MaskPixelType, 3 > MaskImageType;

	typedef itk::Image< PixelType, 3 > OutputImageType;

	typedef itk::DiffusionTensor3D< PixelType > DTITensorType;
	typedef vnl_vector_fixed< PixelType, 15 > DKITensorType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileReader< MaskImageType > MaskReaderType;

	typedef vnl_vector< PixelType > VectorType;
	typedef vnl_matrix< PixelType > MatrixType;

	typedef itk::ImageLinearConstIteratorWithIndex< ImageType > ConstIterator4DType;
	typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType > ConstIterator3DType;
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > Iterator3DType;

	typedef itk::ImageFileWriter< OutputImageType > WriterType;

	typedef boost::tokenizer< boost::char_separator< char > > TokType;
	typedef boost::mt19937 random_number_type;


	/**
	 * Return 2-rand DTI tensor.
	 *
	 * ITK tensor layout =>
	 * 					| 0  1  2  |
	 *			        | X  3  4  |
	 *       			| X  X  5  |
	 */
	static DTITensorType GetDTITensor( const VectorType& theta )
	{
		DTITensorType D( 0.0 );

		D( 0, 0 ) = theta( 1 );
		D( 0, 1 ) = theta( 4 );
		D( 0, 2 ) = theta( 6 );
		D( 1, 1 ) = theta( 2 );
		D( 1, 2 ) = theta( 5 );
		D( 2, 2 ) = theta( 3 );

		return D;
	}

	/**
	 * Return 4-rand DKI tensor (for now a fixed vector with 15 elements).
	 *
	 * theta_dki = [W_1111, W_2222, W_3333,
	 * 				W_1112, W_2223, W_2221,
	 * 				W_2223, W_3331, W_3332,
	 * 				W_1122, W_1133, W_2233,
	 * 				W_1123, W_1223, W_1233]
	 */
	DKITensorType GetDKITensor( const VectorType& theta )
	{
		DKITensorType tensor;
		for( unsigned int i = 0; i < 15; i++ )
			tensor( i ) = theta( 7 + i );

		return tensor;
	}

	/**
	 * Predicate to remove non positive semi-definite DTI/DKI tensors.
	 */
	class non_positive_definite
	{
		public:

		/**
		 * Check for non positive semi-definite tensor.
		 */
		bool operator() ( const VectorType& beta, bool dki_tensor = false )
		{
			DTITensorType DT = GetDTITensor( beta );

			DTITensorType::EigenValuesArrayType eigenValues;
			DT.ComputeEigenValues( eigenValues );

			bool non_negative = false;

			for( unsigned int i = 0; i < 3; i++ )
				if( eigenValues[ i ] < 0 )
					non_negative = true;

			if( dki_tensor )
			{
				DKITensorType DK = GetDKITensor( beta );
				// TODO calculate if DK tensor is postive semi-definite ...
			}

			return non_negative;
		}
	};


	/**
	 * Date, 06-01-2011.
	 *
	 * Implement Gibbs sampler for linear DTI fits.
	 *
	 * Ref: "Introduction to Applied Bayesian Statistics and Estimation for Social Scientists, Scott M. Lynch, 2007. Springer"
	 */
	class DKIFit
    {
        public:

		ImageType::Pointer Input;
		MaskImageType::Pointer Mask;

		MatrixType X_dti; // dti design matrix
		MatrixType X_dki; // dki design matrix
		MatrixType Z; // covariance matrix
		VectorType B; // bvals
		MatrixType R; // bvecs

		OutputImageType::Pointer S0_dti;
		OutputImageType::Pointer S0_dki;

		OutputImageType::Pointer FA_dti;
		OutputImageType::Pointer FA_dki;

		OutputImageType::Pointer Trace_dti;
		OutputImageType::Pointer Trace_dki;

		OutputImageType::Pointer AKC;
		OutputImageType::Pointer Kmax;
		OutputImageType::Pointer Kmin;

		/**
		 * Start DKI fit.
		 */
		void Run( const parameters& args )
        {
			SetImage( args.inputFileName );
			SetMask( args.maskFileName, Input.GetPointer() );
			SetBVals( args.bvalsFileName, args.sep );
			SetBVecs( args.bvecsFileName, args.sep );
			SetDTIDesignMatrix();
			SetDKIDesignMatrix();
			AllocateOutput( Input.GetPointer() );
			Fit( args.iterations, args.regularize );
			Write( args.outputFileName );
        }

        protected:

		/**
		 * Return cholesky decomposition for symmetric matrix.
		 *
		 * Sqrt( A ) using cholesky decomposition is twice as fast as svd/qr.
		 *
		 * The cholesky decomposition decomposes symmetric A = LL'
		 * where L is lower triangular.
		 */
		MatrixType Cholesky( const MatrixType& A )
		{
			/* Test Cholesky ...
			MatrixType A( 3, 3 );
			A( 0, 0 ) = 2; A( 0, 1 ) = 0; A( 0, 2 ) = 0.1;
			A( 1, 0 ) = 0; A( 1, 1 ) = 2; A( 1, 2 ) = 0;
			A( 2, 0 ) = 0.1; A( 2, 1 ) = 0; A( 2, 2 ) = 4;
			MatrixType L = Cholesky( A );
			std::cout << "Cholesky: " << L * L.transpose() << std::endl;
			*/
			vnl_cholesky::Operation op = vnl_cholesky::quiet;

			vnl_cholesky chol( A, op );

			return chol.upper_triangle(); // A = U'U
		}

		/**
		 * Normal distribution sampler.
		 */
		VectorType SampleNormal( random_number_type& ran, unsigned int total, double mean, double sigma )
		{
			using namespace boost;

			// select Gaussian probability distribution
			normal_distribution< double > norm_dist( mean, sigma );

			// bind random number generator to distribution, forming a function
			variate_generator< random_number_type&, normal_distribution< double > > sampler( ran, norm_dist );

			VectorType samples( total );

			for( unsigned int i = 0; i < total; i++ )
				samples( i ) = sampler(); // sample from the distribution

			return samples;
		}

		/**
		 * Inverse gamma distribution sampler.
		 *
		 * Implemented as 1 / gamma, with boost library 1.45 switch to inverse_gamma_distribution()
		 *
		 * shape =~ alpha
		 * scale =~ beta
		 */
		VectorType SampleIGamma( random_number_type& ran, unsigned int total, double alpha, double beta )
		{
			using namespace boost;

		    gamma_distribution< double > gamma_dist( alpha );

		    variate_generator< random_number_type&, gamma_distribution< double > > sampler( ran, gamma_dist );

		    VectorType samples( total );

		    for( unsigned int i = 0; i < total; i++ )
		    	samples( i ) = 1 / ( beta * sampler() ); // sample gamma distribution with scale = 1, and invert

		    return samples;
		}


		/**
		 * Design matrix with ln(S0) as constant .
		 *
		 * See: "C.G. Koay et al / Journal of Magnetic Resonance 128 (2006) pag. 116"
		 *
		 * DTI tensor thus: ln(S0), Dxx, Dyy, Dzz, Dxy, Dyz, Dxz.
		 */
		void SetDTIDesignMatrix()
		{
			MatrixType A( R.rows(), 7, 0 ); // 5 x 7

			for ( unsigned int i = 0; i < R.rows(); i++ )
			{
				double b = B( i );
				double Gx = R( i, 0 );
				double Gy = R( i, 1 );
				double Gz = R( i, 2 );

				A( i, 0 ) = 1; // ln(S0 )

				A( i, 1 ) = -b * Gx * Gx;

				A( i, 2 ) = -b * Gy * Gy;

				A( i, 3 ) = -b * Gz * Gz;

				A( i, 4 ) = -2 * b * Gx * Gy;

				A( i, 5 ) = -2 * b * Gy * Gz;

				A( i, 6 ) = -2 * b * Gx * Gz;
			}

			X_dti = A;
		}

		/**
		 * Bayesian Linear least squares fit.
         *
		 * // TODO add covariance matrix W => to get linear 'weighted' regression.
		 * beta = (X' * X)^-1 * X' * ln( y )
		 */
		VectorType FitVoxel( const MatrixType& X, const VectorType& y, unsigned int iterations,
				unsigned int k, unsigned int n, bool regularize )
		{
			// (X' * X)^-1

			vnl_svd< PixelType > svd( X.transpose() * X );
			MatrixType XTXI = svd.inverse();

			// ln( y )

			VectorType Y( y.size() );
			for( unsigned int i = 0; i < y.size(); i++ )
				Y( i ) = log( y( i ) );

			// (X' * X)^-1 * X' * ln( y )

			VectorType beta = XTXI * X.transpose() * Y;

			// sample sigma from it's inverse gamma marginal distribution

			// Create a Mersenne twister random number generator
			random_number_type eng( static_cast< unsigned int > ( std::time( 0 ) ) );

			PixelType shape = 0.5 * ( n - k );
			PixelType scale = 0.5 * dot_product( ( Y - X * beta ), ( Y - X * beta ) );
			VectorType s2 = SampleIGamma( eng, iterations, shape, scale );

			// -----------------------------------
			// sample distribution for beta (MVN):
			// -----------------------------------
			// mean => inv(X'X) * (X'Y)				[fitted beta]
			// variance => sigma^2 * inv(X'X)		[sqrt using cholesky decomp.]

			std::list< VectorType > bs;
			for( unsigned int i = 0; i < iterations; i++ )
				bs.push_back( beta + SampleNormal( eng, k, 0, 1 ) * Cholesky( s2( i ) * XTXI ) );

			// remove non positive semi definite tensors
			if (regularize )
				bs.remove_if( non_positive_definite() );

			// if nothing left, return zero tensor ...

			if( bs.empty() )
			{
				//std::cout << "Voxel found with sampled tensors not positive semi-definite!" << std::endl;
				return VectorType( beta.size(), 0 ); // return empty tensor
			}

			// get mean tensor

			VectorType mean_beta( beta.size(), 0 );

			for( std::list< VectorType >::iterator it = bs.begin(); it != bs.end(); ++it )
				mean_beta += *it;

			return mean_beta /= bs.size();
		}

		/**
		 * Read bvals from file.
		 */
		void SetBVals( const std::string& bvalsFileName, const std::string& sep )
		{
			// init
			VectorType b( GetNumberOfRows( bvalsFileName ), 0 );

			// open data ...
			std::ifstream in( bvalsFileName.c_str() );
			if ( in.fail() )
			{
				std::cerr << "*** ERROR ***: could not read from: " << bvalsFileName << std::endl;
				exit( EXIT_FAILURE );
			}

			std::string line;
			unsigned int rowIndex = 0;

			while ( getline( in, line ) )
			{
				TokType tok( line, boost::char_separator< char >( sep.c_str() ) );
				try
				{
					for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
					{
						b( rowIndex ) = boost::lexical_cast< PixelType >( *id );
					}
					rowIndex++;
				}
				catch ( boost::bad_lexical_cast& e )
				{
					std::cout << "*** WARNING ***: could not parse " << bvalsFileName << std::endl;
					e.what();
					exit( EXIT_FAILURE );
				}
			}

			in.close();

			B = b;
		}

		void SetBVecs( const std::string& bvecsFileName, const std::string& sep )
		{
			// init
			MatrixType r( GetNumberOfRows( bvecsFileName ), 3, 0 );

			// open data ...
			std::ifstream in( bvecsFileName.c_str() );
			if ( in.fail() )
			{
				std::cerr << "*** ERROR ***: could not read from: " << bvecsFileName << std::endl;
				exit( EXIT_FAILURE );
			}

			std::string line;
			unsigned int rowIndex = 0;

			while ( getline( in, line ) )
			{
				TokType tok( line, boost::char_separator< char >( sep.c_str() ) );
				try
				{
					unsigned int colIndex = 0;
					for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
					{
						r( rowIndex, colIndex ) = boost::lexical_cast< PixelType >( *id );
						colIndex++;
					}
					rowIndex++;
				}
				catch ( boost::bad_lexical_cast& e )
				{
					std::cout << "*** WARNING ***: could not parse " << bvecsFileName << std::endl;
					std::cout << e.what() << std::endl;
					exit( EXIT_FAILURE );
				}
			}

			in.close();

			R = r;
		}

		/**
		 * Build X_dki from bvecs and bvals input files.
		 *
		 * DTI tensor elements:
		 * --------------------
		 * [D_11, D_12, D_13, D_22, D_23, D_33]
		 *
		 * DKI tensor elements:
		 * --------------------
		 * [W_1111, W_2222, W_3333,
		 * 	W_1112, W_2223, W_2221, W_2223, W_3331, W_3332,
		 *  W_1122, W_1133, W_2233,
		 * 	W_1123, W_1223, W_1233]
		 */
		void SetDKIDesignMatrix()
		{
			// check ...
			if( B.size() != R.rows() )
			{
				std::cerr << "*** ERROR ***: bvals and bvecs size does not match!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// init ( e.g. [52 x 22] ) ...
			X_dki = MatrixType( R.rows(), 22, 0 );

			for ( unsigned int i = 0; i < R.rows(); i++ )
			{
				double b = B( i );
				double Gx = R( i, 0 );
				double Gy = R( i, 1 );
				double Gz = R( i, 2 );

				// ln( S0 )

				X_dki( i, 0 ) = 1;

				// DTI tensor

				X_dki( i, 1 ) = -b * Gx * Gx;
				X_dki( i, 2 ) = -b * Gy * Gy;
				X_dki( i, 3 ) = -b * Gz * Gz;
				X_dki( i, 4 ) = -2 * b * Gx * Gy;
				X_dki( i, 5 ) = -2 * b * Gy * Gz;
				X_dki( i, 6 ) = -2 * b * Gx * Gz;

				// DKI tensor

				PixelType c = ( b * b ) / 6;

				X_dki( i, 7 ) = c * Gx * Gx * Gx * Gx; // W_1111
				X_dki( i, 8 ) = c * Gy * Gy * Gy * Gy; // W_2222
				X_dki( i, 9 ) = c * Gz * Gz * Gz * Gz; // W_3333

				X_dki( i, 10 ) = 4 * c * Gx * Gx * Gx * Gy; // W_1112
				X_dki( i, 11 ) = 4 * c * Gx * Gx * Gx * Gy; // W_2223
				X_dki( i, 12 ) = 4 * c * Gx * Gy * Gy * Gy; // W_2221
				X_dki( i, 13 ) = 4 * c * Gy * Gy * Gy * Gz; // W_2223
				X_dki( i, 14 ) = 4 * c * Gx * Gz * Gz * Gz; // W_3331
				X_dki( i, 15 ) = 4 * c * Gy * Gz * Gz * Gz; // W_3332

				X_dki( i, 16 ) = 6 * c * Gx * Gx * Gy * Gy; // W_1122
				X_dki( i, 17 ) = 6 * c * Gx * Gx * Gz * Gz; // W_1133
				X_dki( i, 18 ) = 6 * c * Gy * Gy * Gz * Gz; // W_2233

				X_dki( i, 19 ) = 12 * c * Gx * Gx * Gy * Gz; // W_1123
				X_dki( i, 20 ) = 12 * c * Gx * Gy * Gy * Gz; // W_1223
				X_dki( i, 21 ) = 12 * c * Gx * Gy * Gz * Gz; // W_1233
			}
		}

		/**
		 * Return S0.
		 */
		PixelType GetS0( const VectorType& theta )
		{
			return std::exp( theta( 0 ) );
		}

		/**
		 * Return largest AKC value.
		 *
		 * See:
		 * "Principal invariant and inherent parameters of diffusion kurtosis tensors"
		 * "L. Qi et al / J. Math. Anal. Appl. 349 (2009) 165 - 180."
		 *
		 * return [mean(k), min(k), max(k)]
		 */
		VectorType GetKurtosisValues( const PixelType MD, const VectorType& D_app, const DKITensorType& W )
		{
			// get number of b-zeros ...
			unsigned int b0s = 0;
			for( unsigned int i = 0; i < B.size(); i++ )
				if( B( i ) == 0 )
					b0s++;

			VectorType K( R.rows() - b0s, 0 );

			unsigned int dwi_index = 0;

			for ( unsigned int i = 0; i < R.rows(); i++ )
			{
				if( B( i ) != 0 ) // only dwi ...
				{
					double Gx = R( i, 0 );
					double Gy = R( i, 1 );
					double Gz = R( i, 2 );

					PixelType tmp0 = MD * MD / ( D_app( i ) * D_app( i ) ); // D_app contains b0s ...

					PixelType tmp1 = W( 0 ) * Gx * Gx * Gx * Gx + W( 1 ) * Gy * Gy * Gy * Gy + W( 2 ) * Gz * Gz * Gz * Gz;

					PixelType tmp2 = 4 * ( W( 3 ) * Gx * Gx * Gx * Gy + W( 4 ) * Gx * Gx * Gx * Gy + W( 5 ) * Gx * Gy * Gy * Gy + W( 6 ) * Gy * Gy * Gy * Gz + W( 7 ) * Gx * Gz * Gz * Gz + W( 8 ) * Gy * Gz * Gz * Gz );

					PixelType tmp3 = 6 * ( W( 9 ) * Gx * Gx * Gy * Gy + W( 10 ) * Gx * Gx * Gz * Gz + W( 11 ) * Gy * Gy * Gz * Gz );

					PixelType tmp4 = 12 * ( W( 12 ) * Gx * Gx * Gy * Gz + W( 13 ) * Gx * Gy * Gy * Gz + W( 14 ) * Gx * Gy * Gz * Gz );

					K( dwi_index ) = tmp0 * ( tmp1 + tmp2 + tmp3 + tmp4 );

					dwi_index++;
				}
			}

			VectorType output( 3, 0 ); // average K, min K, max K

			output( 0 ) = K.mean();
			output( 1 ) = K.min_value();
			output( 2 ) = K.max_value();

			return output;
		}

		/**
		 * Return number of training data points.
		 */
		unsigned int GetNumberOfRows( const std::string& inputFile )
		{
			std::ifstream in( inputFile.c_str() );

			if ( in.fail() )
			{
				std::cerr << "*** ERROR ***: could not read from: " << inputFile << std::endl;
				exit( EXIT_FAILURE );
			}
			std::string line;
			unsigned int result = 0;

			while ( getline( in, line ) )
				result++;

			in.close();
			return result;
		}

		/**
		 * MCMC DTI fitting.
		 */
		void Fit( unsigned int iterations, bool regularize )
		{
			ConstIterator4DType it( Input, Input->GetLargestPossibleRegion() );
			ConstIterator3DType mit( Mask, Mask->GetLargestPossibleRegion() );

			Iterator3DType itS0_dti( S0_dti, S0_dti->GetLargestPossibleRegion() );
			Iterator3DType itS0_dki( S0_dki, S0_dki->GetLargestPossibleRegion() );

			Iterator3DType itFA_dti( FA_dti, FA_dti->GetLargestPossibleRegion() );
			Iterator3DType itFA_dki( FA_dki, FA_dki->GetLargestPossibleRegion() );

			Iterator3DType itTrace_dti( Trace_dti, Trace_dti->GetLargestPossibleRegion() );
			Iterator3DType itTrace_dki( Trace_dki, Trace_dki->GetLargestPossibleRegion() );

			Iterator3DType itAKC( AKC, AKC->GetLargestPossibleRegion() );
			Iterator3DType itKmin( Kmin, Kmin->GetLargestPossibleRegion() );
			Iterator3DType itKmax( Kmax, Kmax->GetLargestPossibleRegion() );

			it.SetDirection( 3 );
			it.GoToBegin();
			mit.GoToBegin();

			itS0_dti.GoToBegin();
			itS0_dki.GoToBegin();

			itFA_dti.GoToBegin();
			itFA_dki.GoToBegin();

			itTrace_dti.GoToBegin();
			itTrace_dki.GoToBegin();

			itAKC.GoToBegin();
			itKmin.GoToBegin();
			itKmax.GoToBegin();

			unsigned int totalDWI = ( Input->GetLargestPossibleRegion().GetSize() )[3];

			unsigned int sliceIndex = 0;

			// for each series of DWI ...
			while ( !it.IsAtEnd(), !mit.IsAtEnd() )
			{
				if ( mit.Get() != 0 )
				{
					// slice index
					if( sliceIndex != mit.GetIndex()[ 2 ] )
					{
						std::cout << "Processing slice: " << sliceIndex << std::endl;
						sliceIndex = mit.GetIndex()[ 2 ];
					}

					VectorType y( totalDWI );

					while ( !it.IsAtEndOfLine() )
					{
						y( ( it.GetIndex() )[3] ) = it.Get();
						++it;
					}

					// DTI => 7 elements (ln(S0) = 1, theta_D = 6)
					// DTI fit only ...

					unsigned int k = 7;
					unsigned int n = y.size();
					VectorType theta_dti = FitVoxel( X_dti, y, iterations, k, n, regularize );
					itS0_dti.Set( GetS0( theta_dti ) );
					DTITensorType tensor_simple = GetDTITensor( theta_dti );
					itFA_dti.Set( tensor_simple.GetFractionalAnisotropy() );
					itTrace_dti.Set( tensor_simple.GetTrace() );

					// DKI => 22 elements (ln(S0) = 1, theta_D = 6, theta_K = 15)
					// DTI and DKI fit in one...

					k += 15;

					VectorType theta_dki = FitVoxel( X_dki, y, iterations, k, n, regularize );
					itS0_dki.Set( GetS0( theta_dki ) );

					DTITensorType tensor_dti = GetDTITensor( theta_dki ); // DTI tensor from DKI fit ...
					itFA_dki.Set( tensor_dti.GetFractionalAnisotropy() );
					itTrace_dki.Set( tensor_dti.GetTrace() );

					DKITensorType tensor_dki = GetDKITensor( theta_dki );
					VectorType tmp = GetKurtosisValues( tensor_dti.GetTrace(), FitDapp( y ), tensor_dki );

					itAKC.Set( tmp( 0 ) );
					itKmin.Set( tmp( 1 ) );
					itKmax.Set( tmp( 2 ) );
				}

				it.NextLine();
				++mit;

				++itS0_dti;
				++itS0_dki;

				++itFA_dti;
				++itFA_dki;

				++itTrace_dti;
				++itTrace_dki;

				++itAKC;
				++itKmin;
				++itKmax;
			}
		}

		/**
		 * Return apparant diffusion per direction.
		 *
		 * ln( Sb / S0 ) / -b = D_app
		 */
		VectorType FitDapp( const VectorType& y )
		{
			// average S0 ...

			PixelType S0 = 0;
			PixelType totalS0;

			for( unsigned int i = 0; i < B.size(); i++ )
			{
				if( B( i ) == 0 )
				{
					S0 += y( i );
					totalS0++;
				}
			}
			if( totalS0 == 0 )
			{
				std::cerr << "*** ERROR ***: no b-zeros!" << std::endl;
				exit( EXIT_FAILURE );
			}
			S0 /= totalS0;

			// D_app ...

			VectorType D_app( R.rows(), 0 );

			for( unsigned int i = 0; i < B.size(); i++ )
			{
				if( B( i ) != 0 )
					D_app( i ) = log( y( i )  / S0 ) / - B( i );
			}

			return D_app;
		}

		/**
		 * Return covariate matrix Z.
		 */
		MatrixType BuildCovarianceMatrix( const VectorType& y )
		{
			MatrixType Z( y.size(), y.size(), 0 );

			// fill diagonal
			for( unsigned int i = 0; i < Z.rows(); i++ )
				Z( i, i ) = y( i ) * y( i );

			return Z;
		}

		/**
		 * Allocate all output images to 0.
		 */
		void AllocateOutput( ImageType::ConstPointer input )
		{
			MaskImageType::Pointer output = MaskImageType::New();

			OutputImageType::RegionType region3D;
			OutputImageType::IndexType index3D;
			OutputImageType::SizeType size3D;
			OutputImageType::SpacingType spacing3D;
			OutputImageType::PointType origin3D;

			ImageType::RegionType region = input->GetLargestPossibleRegion();
			ImageType::SizeType size = region.GetSize();
			ImageType::IndexType index = region.GetIndex();
			ImageType::SpacingType spacing = input->GetSpacing();
			ImageType::PointType origin = input->GetOrigin();

			size3D[0] = size[0];
			size3D[1] = size[1];
			size3D[2] = size[2];
			index3D[0] = index[0];
			index3D[1] = index[1];
			index3D[2] = index[2];
			origin3D[0] = origin[0];
			origin3D[1] = origin[1];
			origin3D[2] = origin[2];
			spacing3D[0] = spacing[0];
			spacing3D[1] = spacing[1];
			spacing3D[2] = spacing[2];

			region3D.SetSize( size3D );
			region3D.SetIndex( index3D );

			// create
			S0_dti = OutputImageType::New();
			S0_dki = OutputImageType::New();

			FA_dti = OutputImageType::New();
			FA_dki = OutputImageType::New();

			Trace_dti = OutputImageType::New();
			Trace_dki = OutputImageType::New();

			AKC = OutputImageType::New();
			Kmax = OutputImageType::New();
			Kmin = OutputImageType::New();

			// set
			S0_dti->SetRegions( region3D );
			S0_dti->SetSpacing( spacing3D );
			S0_dti->SetOrigin( origin3D );
			S0_dti->Allocate();
			S0_dti->FillBuffer( 0 );

			S0_dki->SetRegions( region3D );
			S0_dki->SetSpacing( spacing3D );
			S0_dki->SetOrigin( origin3D );
			S0_dki->Allocate();
			S0_dki->FillBuffer( 0 );

			FA_dti->SetRegions( region3D );
			FA_dti->SetSpacing( spacing3D );
			FA_dti->SetOrigin( origin3D );
			FA_dti->Allocate();
			FA_dti->FillBuffer( 0 );

			FA_dki->SetRegions( region3D );
			FA_dki->SetSpacing( spacing3D );
			FA_dki->SetOrigin( origin3D );
			FA_dki->Allocate();
			FA_dki->FillBuffer( 0 );

			Trace_dti->SetRegions( region3D );
			Trace_dti->SetSpacing( spacing3D );
			Trace_dti->SetOrigin( origin3D );
			Trace_dti->Allocate();
			Trace_dti->FillBuffer( 0 );

			Trace_dki->SetRegions( region3D );
			Trace_dki->SetSpacing( spacing3D );
			Trace_dki->SetOrigin( origin3D );
			Trace_dki->Allocate();
			Trace_dki->FillBuffer( 0 );

			AKC->SetRegions( region3D );
			AKC->SetSpacing( spacing3D );
			AKC->SetOrigin( origin3D );
			AKC->Allocate();
			AKC->FillBuffer( 0 );

			Kmax->SetRegions( region3D );
			Kmax->SetSpacing( spacing3D );
			Kmax->SetOrigin( origin3D );
			Kmax->Allocate();
			Kmax->FillBuffer( 0 );

			Kmin->SetRegions( region3D );
			Kmin->SetSpacing( spacing3D );
			Kmin->SetOrigin( origin3D );
			Kmin->Allocate();
			Kmin->FillBuffer( 0 );
		}

		/**
		 * Set input image.
		 */
		void SetImage( const std::string& inputFileName )
		{
			if ( !inputFileName.empty() )
			{
				ReaderType::Pointer reader = ReaderType::New();
				reader->SetFileName( inputFileName );
				reader->Update();
				Input = reader->GetOutput();
			}
			else
			{
				std::cerr << "Could not read input!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		/**
		 * Return mask if file given, else create empty mask from input file.
		 */
		void SetMask( const std::string& maskFileName, ImageType::ConstPointer input )
		{
			if ( !maskFileName.empty() )
			{
				MaskReaderType::Pointer reader = MaskReaderType::New();
				reader->SetFileName( maskFileName );
				reader->Update();
				Mask = reader->GetOutput();
			}
			else
			{
				MaskImageType::Pointer output = MaskImageType::New();
				MaskImageType::RegionType region3D;
				MaskImageType::IndexType index3D;
				MaskImageType::SizeType size3D;
				MaskImageType::SpacingType spacing3D;
				MaskImageType::PointType origin3D;

				ImageType::RegionType region = input->GetLargestPossibleRegion();
				ImageType::SizeType size = region.GetSize();
				ImageType::IndexType index = region.GetIndex();
				ImageType::SpacingType spacing = input->GetSpacing();
				ImageType::PointType origin = input->GetOrigin();

				size3D[0] = size[0];
				size3D[1] = size[1];
				size3D[2] = size[2];
				index3D[0] = index[0];
				index3D[1] = index[1];
				index3D[2] = index[2];
				origin3D[0] = origin[0];
				origin3D[1] = origin[1];
				origin3D[2] = origin[2];
				spacing3D[0] = spacing[0];
				spacing3D[1] = spacing[1];
				spacing3D[2] = spacing[2];

				region3D.SetSize( size3D );
				region3D.SetIndex( index3D );

				// set
				output->SetRegions( region3D );
				output->SetSpacing( spacing3D );
				output->SetOrigin( origin3D );
				output->Allocate();
				output->FillBuffer( 1 );

				Mask = output;
			}
		}

		/**
		 * Write output images.
		 */
		void Write( const std::string& outputFileName )
		{
			// Writer
			WriterType::Pointer writer = WriterType::New();

			// S0
			writer->SetFileName( ( outputFileName + std::string( "_S0_dti.nii.gz" ) ).c_str() );
			writer->SetInput( S0_dti );
			writer->Update();

			writer->SetFileName( ( outputFileName + std::string( "_S0_dki.nii.gz" ) ).c_str() );
			writer->SetInput( S0_dki );
			writer->Update();

			// FA
			writer->SetFileName( ( outputFileName + std::string( "_FA_dti.nii.gz" ) ).c_str() );
			writer->SetInput( FA_dti );
			writer->Update();

			writer->SetFileName( ( outputFileName + std::string( "_FA_dki.nii.gz" ) ).c_str() );
			writer->SetInput( FA_dki );
			writer->Update();

			// Trace
			writer->SetFileName( ( outputFileName + std::string( "_Trace_dti.nii.gz" ) ).c_str() );
			writer->SetInput( Trace_dti );
			writer->Update();

			writer->SetFileName( ( outputFileName + std::string( "_Trace_dki.nii.gz" ) ).c_str() );
			writer->SetInput( Trace_dki );
			writer->Update();

			// AKC
			writer->SetFileName( ( outputFileName + std::string( "_AKC.nii.gz" ) ).c_str() );
			writer->SetInput( AKC );
			writer->Update();

			// Kmin
			writer->SetFileName( ( outputFileName + std::string( "_Kmin.nii.gz" ) ).c_str() );
			writer->SetInput( Kmin );
			writer->Update();

			// Kmax
			writer->SetFileName( ( outputFileName + std::string( "_Kmax.nii.gz" ) ).c_str() );
			writer->SetInput( Kmax );
			writer->Update();
		}

    }; // end class DKIFit

} // end namespace dki


/**
 * DKI fitting routine.
 */
int main( int argc, char ** argv )
{
    tkd::CmdParser p( argv[0], "Fit voxel-wise DTI using Gibbs sampling and postitive semi-definite tensor regularization." );

    dki::parameters args;
    args.sep = " ";
    args.iterations = 1000;
    args.regularize = false;

    p.AddArgument( args.inputFileName, "input" )
        ->AddAlias( "i" )
        ->SetDescription( "Input 4D image" )
        ->SetRequired( true );

    p.AddArgument( args.bvecsFileName, "bvecs" )
        ->AddAlias( "r" )
        ->SetDescription( "Diffusion gradient vectors (FSL format)" )
        ->SetRequired( true );

    p.AddArgument( args.bvalsFileName, "bvals" )
        ->AddAlias( "b" )
        ->SetDescription( "Diffusion b-value vector (FSL format)" )
        ->SetRequired( true );

    p.AddArgument( args.outputFileName, "output" )
        ->AddAlias( "o" )
        ->SetDescription( "Output filename base" )
        ->SetRequired( true );

    p.AddArgument( args.maskFileName, "mask" )
        ->AddAlias( "m" )
        ->SetDescription( "Mask 3D image" );

    p.AddArgument( args.sep, "separation" )
        ->AddAlias( "s" )
        ->SetDescription( "Separation string in bvecs file (default: ' ')" );

    p.AddArgument( args.iterations, "iterations" )
        ->AddAlias( "it" )
        ->SetDescription( "Number of MCMC sampling iterations per voxel (default: 1000)" );

    p.AddArgument( args.regularize, "regularize" )
        ->AddAlias( "reg" )
        ->SetDescription( "Regularize posteriors by removal of non postivive semi-definite tensors (default: false)" );

    if ( !p.Parse( argc, argv ) )
    {
        p.PrintUsage( std::cout );
        return EXIT_FAILURE;
    }

    dki::DKIFit fit;
    fit.Run( args );

    return EXIT_SUCCESS;
}





