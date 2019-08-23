#include "tkdCmdParser.h"

#include "vnl/algo/vnl_svd.h"
#include "vnl_vector_to_std_vector.h"

#include <algorithm>
#include "dkifitCommon.h"

#include "quantile.hpp"




namespace dki
{

	/**
	 * ------------------------------------------------------------------------------------
	 * Fit linear diffusion kurtosis imaging tensors using RESTORE.
	 *
	 * Robust Diffusion Tensor Estimation.
	 *
	 * Magnetic Resonance in Medicine 53:1088-1095 (2005)
	 * Lin-Ching Chang, Derek K. Jones and Carlo Pierpaoli
	 * ------------------------------------------------------------------------------------
	 */
	class RESTOREDkiFit : public DkiFit
	{

	public:

		const static double THRESHOLD_CONVERGENCE = 0.0001;
		const static double MULTIPLY_NOISE = 3;

		/**
		 * Extended arguments.
		 */
		struct restore_parameters : parameters
		{
			PixelType sigma;
			int max_iter;
		};

		/**
		 * Extended constructor.
		 */
		RESTOREDkiFit( const restore_parameters& args ) : DkiFit( args )
		{
			m_max_iter = args.max_iter;
			m_sigma = args.sigma;

			FitAllVoxels( args );
			Write( args.outputFileName );
		}

	protected:

		unsigned int m_max_iter;
		PixelType m_sigma;

		/**
		 * Fit for every voxel in mask.
		 */
		void FitAllVoxels( const parameters& args )
		{
			ConstIterator4DType it( m_Input, m_Input->GetLargestPossibleRegion() );
			Iterator4DType ot( m_Output, m_Output->GetLargestPossibleRegion() );
			ConstIterator3DType mit( m_Mask, m_Mask->GetLargestPossibleRegion() );

			it.SetDirection( 3 );
			it.GoToBegin();

			ot.SetDirection( 3 );
			ot.GoToBegin();

			mit.GoToBegin();

			unsigned int sliceIndex = 0;

			unsigned int totalDWI = ( m_Input->GetLargestPossibleRegion().GetSize() )[3];

			// for each series of DWI ...
			while ( !it.IsAtEnd(), !mit.IsAtEnd() )
			{
				if ( mit.Get() != 0 )
				{
					// slice index
					if ( sliceIndex != mit.GetIndex()[2] )
					{
						std::cout << "Processing slice: " << sliceIndex << std::endl;
						sliceIndex = mit.GetIndex()[2];
					}

					VectorType allVoxelData( totalDWI );

					while ( !it.IsAtEndOfLine() )
					{
						allVoxelData( ( it.GetIndex() )[3] ) = it.Get();
						++it;
					}

					 // data = ln( S_M / S_0 )

					VectorType b = LogDWIOnlyData( allVoxelData );

					VectorType x = RESTORE( b );

					for ( unsigned int i = 0; i < x.size(); i++ )
					{
						ot.Set( x( i ) );
						++ot;
					}
				}

				it.NextLine();
				ot.NextLine();
				++mit;
			}
		}

		/**
		 * Return variance matrix with 1 / w^2 on diagonal.
		 */
		DiagMatrixType VarianceMatrixConstantWeights( unsigned int size, PixelType w )
		{
			DiagMatrixType W( size, size );
			for( unsigned int i = 0; i < size; i++ )
				W( i, i ) = 1. / w * w;

			return W;
		}

		/**
		 * Unconstrained weighted linear least squares using RESTORE.
		 */
		VectorType RESTORE( const VectorType& b )
		{
			// multiply background SD: Henkelman RM. Med Phys 1985; 12:232-233

			DiagMatrixType W = VarianceMatrixConstantWeights( m_M, m_sigma * 1.5267 );

			vnl_svd< double > inverter( CheckDeterminant( m_DesignMatrix.transpose() * W * m_DesignMatrix ) );
			MatrixType AWAinv = inverter.inverse();
			VectorType x = AWAinv * m_DesignMatrix.transpose() * W * b;

			VectorType r = Residuals( m_DesignMatrix, W, b, x );

			IndexVectorType outliers = Outliers( r );

			if( outliers.sum() > 0 ) // iterative reweighting process
			{
				GMM( W, b, x );

				VectorType r = Residuals( m_DesignMatrix, W, b, x );
				outliers = Outliers( r );

				if( outliers.sum() > 0 )
				{
					return FitWithoutOutliers( outliers, b );
				}
			}

			return x;
		}

		/**
		 * Remove outlier entries from A and b
		 */
		VectorType FitWithoutOutliers( const IndexVectorType& outliers, const VectorType& b )
		{
			unsigned int rows = outliers.size() - outliers.sum();

			MatrixType reducedA = MatrixType( rows, m_DesignMatrix.cols() );
			VectorType reducedb = VectorType( rows );

			unsigned int ind = 0;

			for( unsigned int i = 0; i < outliers.size(); i++ )
				if( outliers( i ) == 0 )
				{
					reducedA.set_row( ind, m_DesignMatrix.get_row( i ) );
					reducedb( ind ) = b( i );
					ind++;
				}

			DiagMatrixType W = VarianceMatrixConstantWeights( m_M, m_sigma * 1.5267 );
			vnl_svd< double > inverter( CheckDeterminant( reducedA.transpose() * W * reducedA ) );
			MatrixType AWAinv = inverter.inverse();
			VectorType x = AWAinv * reducedA.transpose() * W * reducedb;

			return x;
		}

		/**
		 * Return vector with value 1 if outlier, 0 otherwise.
		 *
		 * x is outlier if:
		 *
		 * 	x < Q1 - 1.5(IQR)
		 * 	x > Q3 + 1.5(IQR)
		 */
		IndexVectorType Outliers( const VectorType& r )
		{
			std::vector< PixelType > data = vnl_vector_to_std_vector< PixelType >( r );

			// sort min - max

			std::sort( data.begin(), data.end() );

			PixelType Q1 = boost::svg::quantile( data, 0.25, 4 );
			PixelType Q3 = boost::svg::quantile( data, 0.75, 4 );

			PixelType IQR = Q3 - Q1;
			PixelType lowerBound = Q1 - 1.5 * IQR;
			PixelType upperBound = Q3 + 1.5 * IQR;

			IndexVectorType outliers( r.size(), 0 );

			for( unsigned int i = 0; i < r.size(); i++ )
				if( ( r( i ) < lowerBound ) || ( r( i ) > upperBound ) )
					outliers( i ) = 1;

			return outliers;
		}

		/**
		 *
		 * German-McClure M-estimator (GMM): iterative reweighted least-squares fitting.
		 *
		 * Update weights until convergence (is insignificant changes in the largest residual).
		 *
		 * Return vector with reweighted x.
		 */
		VectorType GMM( const DiagMatrixType& Wstart, const VectorType& b, const VectorType& xStart )
		{
			unsigned int iter = 0;
			PixelType error = 1.0;

			VectorType xOld = xStart;
			VectorType xNew = xStart;
			DiagMatrixType W = Wstart;

			while ( error > THRESHOLD_CONVERGENCE )
			{
				xOld = xNew;

				if( iter > m_max_iter )
				{
					std::cerr << "*** ERROR ***: maximum iteration reached. No convergence."<< std::endl;
					return( xNew );
				}

				// weighted linear least square fit

				MatrixType AWA = m_DesignMatrix.transpose() * W * m_DesignMatrix;

				vnl_svd< PixelType > svd( CheckDeterminant( AWA ) );

				MatrixType AWAinv = svd.inverse();
				xNew = AWAinv * m_DesignMatrix * W * b;

				// Reweight residuals

				VectorType residuals = Residuals( m_DesignMatrix, W, b, xNew );

				W = GMMWeightFunctionMatrix( residuals );

				error = ( m_DesignMatrix * xOld - m_DesignMatrix * xNew ).max_value();
				iter++;
			}

			return xNew;
		}

		/**
		 * GMM weight function matrix: 1 / r^2 + C^2
		 */
		DiagMatrixType GMMWeightFunctionMatrix( const VectorType& residuals )
		{
			DiagMatrixType W( residuals.size(), 0 );

			PixelType C = GMMScalingFactorC( residuals );

			for( unsigned int i = 0; i < residuals.size(); i++ )
				if( residuals( i ) != 0 )
					W( i, i ) = 1 / residuals( i ) * residuals( i ) * C * C;
				else
					W( i, i ) = 1.0;

			return W;
		}

		/**
		 * Return scaling factor C.
		 */
		PixelType GMMScalingFactorC( const VectorType& V )
		{
			PixelType median = Median( V );

			VectorType diff( V );
			for( unsigned int i = 0; i < V.size(); i++ )
				diff = vcl_abs( V( i ) - median );

			return Median( diff );
		}


	}; // end class RESTOREDkiFit

} // end namespace dki


/**
 * RESTORE dkifit
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "RESTORE diffusion kurtosis tensor estimation" );

	dki::RESTOREDkiFit::restore_parameters args;

	args.sep = " ";
	args.max_iter = 1000;

	p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input 4D image" ) ->SetRequired( true );

	p.AddArgument( args.bvecsFileName, "bvecs" ) ->AddAlias( "r" ) ->SetDescription( "Diffusion gradient vectors (FSL format)" ) ->SetRequired(
			true );

	p.AddArgument( args.bvalsFileName, "bvals" ) ->AddAlias( "b" ) ->SetDescription( "Diffusion b-value vector (FSL format)" ) ->SetRequired(
			true );

	p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output tensor filename" ) ->SetRequired( true );

	p.AddArgument( args.maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask 3D image" );

	p.AddArgument( args.sep, "separation" ) ->AddAlias( "s" ) ->SetDescription( "Separation string in bvecs file (default: ' ')" );

	p.AddArgument( args.sigma, "sigma" ) ->SetDescription( "SD of the background noise" ) ->SetRequired( true );

	p.AddArgument( args.max_iter, "iterations" ) ->AddAlias( "iter" ) ->SetDescription( "Maximum iterations for IRLS (default: 1000)" );



	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	dki::RESTOREDkiFit fit( args );

	return EXIT_SUCCESS;
}
