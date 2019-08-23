#include "tkdCmdParser.h"

#include "vnlquadprog.hpp"

#include "dkifitCommon.h"


namespace dki
{

	/**
	 * Fit constrained diffusion kurtosis imaging tensors by solving the quadratic
	 * programming problem with linear constraint (CLLS-QP).
	 * ------------------------------------------------------------------------------------
	 * "Estimation of Tensors and Tensor-Derived Measures in Diffusional Kurtosis Imaging".
	 * Ali Tabesh, Jens H. Jensen, Babak A. Ardekani and Joseph A. Helpern
	 *
	 * Magnetic Resonance in Medicine 000:000-000 (2010).
	 * ------------------------------------------------------------------------------------
	 * Part of this code was first implemented in matlab (Umesh: fit_dki_lin.m 23-03-2011).
	 */
	class ConstrainedDkiFit : public DkiFit
	{

	public:

		/**
		 * Extended arguments.
		 */
		struct constrained_parameters : parameters
		{
			int max_iter;
			double max_kurtosis;
		};


		PixelType m_max_kurtosis;
		unsigned int m_max_iter;
		MatrixType m_f_partial;
		MatrixType m_Hessian; // A'A / 2, where A = design matrix
		MatrixType m_ConstrainsMatrix;
		VectorType m_d; // zeros


		/**
		 * Extended constructor.
		 */
		ConstrainedDkiFit( const constrained_parameters& args ) : DkiFit( args )
		{

			m_max_kurtosis = args.max_kurtosis;
			m_max_iter = args.max_iter;

			InitHessian_and_fpartial();
			InitConstrainsMatrix();
			Initd();

			FitAllVoxels( args );

			Write( args.outputFileName );
		}

	protected:

		/**
		 * d as zeros in: Cx <= d.
		 */
		void Initd()
		{
			m_d = VectorType( m_ConstrainsMatrix.rows(), 0 );
		}

		/**
		 * Initialize the f_part and Hessian matrix.
		 *
		 * f_part = - 2 * transpose(DesignMatrix);
		 * H = 2 * transpose(DesignMatrix) * DesignMatrix;
		 */
		void InitHessian_and_fpartial()
		{
			m_f_partial = -2. * m_DesignMatrix.transpose();
			m_Hessian = 2. * m_DesignMatrix.transpose() * m_DesignMatrix;
		}

		/**
		 * Initialize the constrains matrix C ((13) in Tabesh et al).
		 *
		 * 						| -A_D				 0 	 |
		 * ConstrainsMatrix C = |  0		 		-A_K |
		 * 						| -(C/b_max)A_D		 A_K |
		 */
		void InitConstrainsMatrix()
		{
			unsigned int rows = m_A_D.rows() + 2 * m_A_K.rows();
			unsigned int cols = m_A_D.cols() + m_A_K.cols();

			m_ConstrainsMatrix = MatrixType( rows, cols, 0 );

			PixelType C = - m_max_kurtosis / m_bvals.max_value();

			// upper corner left (-A_D)
			for( unsigned int r = 0; r < m_A_D.rows(); r++ )
				for( unsigned int c = 0; c < m_A_D.cols(); c++ )
					m_ConstrainsMatrix( r, c ) = - m_A_D( r, c );

			// lower corner left (-(C/b_max)A_D)
			for( unsigned int r = 0; r < m_A_D.rows(); r++ )
				for( unsigned int c = 0; c < m_A_D.cols(); c++ )
					m_ConstrainsMatrix( r + m_A_D.rows() + m_A_K.rows(), c ) = C * m_A_D( r, c );

			// middle corner right (-A_K)
			for( unsigned int r = 0; r < m_A_K.rows(); r++ )
				for( unsigned int c = 0; c < m_A_K.cols(); c++ )
					m_ConstrainsMatrix(r + m_A_D.rows(), c + m_A_D.cols()) = - m_A_K( r, c );

			// lower corner right (A_K)
			for( unsigned int r = 0; r < m_A_K.rows(); r++ )
				for( unsigned int c = 0; c < m_A_K.cols(); c++ )
					m_ConstrainsMatrix(r + m_A_D.rows() + m_A_K.rows(), c + m_A_D.cols()) = m_A_K( r, c );
		}

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

					VectorType param = CLLS_QP( LogDWIOnlyData( allVoxelData ) );

					for ( unsigned int i = 0; i < param.size(); i++ )
					{
						ot.Set( param( i ) );
						++ot;
					}
				}

				it.NextLine();
				ot.NextLine();
				++mit;
			}
		}

		/**
		 * Constrained linear least squares quadratic programming using uBLAS.
		 */
		VectorType CLLS_QP( const VectorType& data )
		{
			VectorType f = m_f_partial * data;
			VectorType x( m_Hessian.rows(), 0 );

			// NO equality constraints (i.e. Ax = b), only inequality constraints(i.e. Cx <= d)

			qp::solve_quadprog( m_Hessian, f, MatrixType(0,0), VectorType(0), m_ConstrainsMatrix, m_d, x );

			return x;
		}



	}; // end class ConstrainedDkiFit

} // end namespace dki


/**
 * Constrained dkifit.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Constrained diffusion kurtosis tensor estimation (CLSS-QP according Tabesh et al MRM 2010)" );

	dki::ConstrainedDkiFit::constrained_parameters args;

	args.sep = " ";
	args.max_kurtosis = 3.0;
	args.max_iter = 200;

	p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input 4D image" ) ->SetRequired( true );

	p.AddArgument( args.bvecsFileName, "bvecs" ) ->AddAlias( "r" ) ->SetDescription( "Diffusion gradient vectors (FSL format)" ) ->SetRequired(
			true );

	p.AddArgument( args.bvalsFileName, "bvals" ) ->AddAlias( "b" ) ->SetDescription( "Diffusion b-value vector (FSL format)" ) ->SetRequired(
			true );

	p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output tensor filename" ) ->SetRequired( true );

	p.AddArgument( args.maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask 3D image" );

	p.AddArgument( args.sep, "separation" ) ->AddAlias( "s" ) ->SetDescription( "Separation string in bvecs file (default: ' ')" );

	p.AddArgument( args.max_kurtosis, "max-kurtosis" ) ->AddAlias( "C" ) ->SetDescription( "Maximum kurtosis value (default: 3)" );

	p.AddArgument( args.max_iter, "max-iterations" ) ->AddAlias( "mit" ) ->SetDescription( "Maximum number of iterations (default: 200)" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	dki::ConstrainedDkiFit fit( args );

	return EXIT_SUCCESS;
}
