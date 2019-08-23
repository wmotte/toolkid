#include "tkdCmdParser.h"

#include "vnl/algo/vnl_svd.h"

#include "dkifitCommon.h"


namespace dki
{

	/**
	 * ------------------------------------------------------------------------------------
	 * Fit linear diffusion kurtosis imaging tensors.
	 * ------------------------------------------------------------------------------------
	 */
	class LinearDkiFit : public DkiFit
	{

	public:

		/**
		 * Extended arguments.
		 */
		struct linear_parameters : parameters
		{
		};

		/**
		 * Extended constructor.
		 */
		LinearDkiFit( const linear_parameters& args ) : DkiFit( args )
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

					VectorType x = ULLS( b );

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
		 * Unconstrained linear least squares.
		 *
		 * x = pinv(A) * b
		 */
		VectorType ULLS( const VectorType& b )
		{
			VectorType x( dki::DkiFit::m_OUTPUT_IMAGES, 0 );

			vnl_svd< double > inverter( m_DesignMatrix );
			x = inverter.solve( b );

			return x;
		}

	}; // end class LinearDkiFit

} // end namespace dki


/**
 * Linear dkifit
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Linear diffusion kurtosis tensor estimation" );

	dki::LinearDkiFit::linear_parameters args;

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

	dki::LinearDkiFit fit( args );

	return EXIT_SUCCESS;
}
