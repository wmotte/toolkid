#include "tkdCmdParser.h"

#include "vnl/algo/vnl_amoeba.h"

#include "dkiOptimizeFunction.h"
#include "dkifitCommon.h"

namespace dki
{
	/**
	 * ------------------------------------------------------------------------------------
	 * Fit nonlinear amoeba diffusion kurtosis imaging tensors.
	 * ------------------------------------------------------------------------------------
	 */
	class NonlinearDkiFit : public DkiFit
	{

	public:

		typedef DKIOptimizeFunction FunctionType;
		typedef vnl_amoeba OptimizerType;

		/**
		 * Extended arguments.
		 */
		struct nonlinear_parameters : parameters
		{
		};

		/**
		 * Extended constructor.
		 */
		NonlinearDkiFit( const nonlinear_parameters& args ) : DkiFit( args )
		{
			FitAllVoxels( args );
			Write( args.outputFileName );
		}

	protected:

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

					VectorType x = UNLS( b );

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
		 * Unconstrained nonlinear least squares.
		 */
		VectorType UNLS( const VectorType& b )
		{
			VectorType x( dki::DkiFit::m_OUTPUT_IMAGES, 0 );

			// initial values

			vnl_svd< double > inverter( m_DesignMatrix );
			x = inverter.solve( b );


			FunctionType function( m_DesignMatrix, b );
			OptimizerType optimizer( function );

			// minimize

			optimizer.minimize( x );

			return x;
		}
	}; // end class NonlinearDkiFit

} // end namespace dki


/**
 * Nonlinear dkifit
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Nonlinear diffusion kurtosis tensor estimation" );

	dki::NonlinearDkiFit::nonlinear_parameters args;

	args.sep = " ";

	p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input 4D image" ) ->SetRequired( true );

	p.AddArgument( args.bvecsFileName, "bvecs" ) ->AddAlias( "r" ) ->SetDescription( "Diffusion gradient vectors (FSL format)" ) ->SetRequired(
			true );

	p.AddArgument( args.bvalsFileName, "bvals" ) ->AddAlias( "b" ) ->SetDescription( "Diffusion b-value vector (FSL format)" ) ->SetRequired(
			true );

	p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output tensor filename" ) ->SetRequired( true );

	p.AddArgument( args.maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask 3D image" );

	p.AddArgument( args.sep, "separation" ) ->AddAlias( "s" ) ->SetDescription( "Separation string in bvecs file (default: ' ')" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	dki::NonlinearDkiFit fit( args );

	return EXIT_SUCCESS;
}

