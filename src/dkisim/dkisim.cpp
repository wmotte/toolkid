#include "tkdCmdParser.h"

#include "dkifitCommon.h"
#include <limits>

namespace dki
{

	/**
	 * Simulate kurtosis data by means of two crossing-fiber populations.
	 *
	 * Init from DkiFit to read bvecs.
	 */
	class DkiSim : public DkiFit
	{

	public:

		/**
		 * Extended arguments.
		 */
		struct sim_parameters : parameters
		{
			std::vector< int > size3D;
			std::vector< PixelType > simBvals;
			int numberOfBzeros;
			int numberOfOutliers;
			PixelType SNR;
			PixelType B0Intensity;
			PixelType FA;
			PixelType MD;
			PixelType el1;
			PixelType el2;
			PixelType az1;
			PixelType az2;
			PixelType fraction1;
			PixelType fraction2;
		};

		/**
		 * Constructor.
		 */
		DkiSim( const sim_parameters& args )
		{
			m_sep = args.sep;
			m_simBvals = args.simBvals;

			m_numberOfBzeros = args.numberOfBzeros;
			m_numberOfOutliers = args.numberOfOutliers;
			m_SNR = args.SNR;
			m_B0Intensity = args.B0Intensity;
			m_FA = args.FA;
			m_MD = args.MD;
			m_el1 = args.el1;
			m_el2 = args.el2;
			m_az1 = args.az1;
			m_az2 = args.az2;
			m_fraction1 = args.fraction1;
			m_fraction2 = args.fraction2;

			ReadBVecs( args.bvecsFileName );
			RemainDWIOnly();
			InitData( args.size3D );
			FillVoxels();
			Write( args.outputFileName );
		}

	protected:

		std::vector< int > m_outputSize;
		std::vector< PixelType > m_simBvals;
		int m_numberOfBzeros;
		int m_numberOfOutliers;
		PixelType m_SNR;
		PixelType m_B0Intensity;

		PixelType m_FA;
		PixelType m_MD;

		PixelType m_el1;
		PixelType m_el2;

		PixelType m_az1;
		PixelType m_az2;

		PixelType m_fraction1;
		PixelType m_fraction2;


		/**
		 * Fit for every voxel in mask.
		 */
		void FillVoxels()
		{
			Iterator4DType it( m_Output, m_Output->GetLargestPossibleRegion() );
			it.GoToBegin();
			it.SetDirection( 3 );

			// for each 4D voxel

			while ( !it.IsAtEnd() )
			{

				VectorType S = SimulateDiffusionTensorSignal();
				VectorType Sn = AddNoiseToDWIs();

				for ( unsigned int i = 0; i < Sn.size(); i++ )
				{
					it.Set( Sn( i ) );
					++it;
				}

				it.NextLine();
			}
		}

		/**
		 * TODO
		 */
		VectorType SimulateDiffusionTensorSignal()
		{
			return VectorType( m_M, 1 );
		}

		/**
		 * TODO
		 */
		VectorType AddNoiseToDWIs()
		{
			return VectorType( m_M, 1.2 );
		}

		/**
		 * Init data (read input, allocate output).
		 */
		void InitData( const std::vector< int >& size3D )
		{
			OutputImageType::RegionType region;
			OutputImageType::SizeType size;
			OutputImageType::IndexType index;
			OutputImageType::SpacingType spacing;
			OutputImageType::PointType origin;

			size[ 0 ] = size3D[ 0 ];
			size[ 1 ] = size3D[ 1 ];
			size[ 2 ] = size3D[ 2 ];
			size[ 3 ] = m_M + m_numberOfBzeros;

			index[ 0 ] = 0;
			index[ 1 ] = 0;
			index[ 2 ] = 0;
			index[ 3 ] = 0;

			spacing[ 0 ] = 1.0;
			spacing[ 1 ] = 1.0;
			spacing[ 2 ] = 1.0;
			spacing[ 3 ] = 1.0;

			origin[ 0 ] = 0;
			origin[ 1 ] = 0;
			origin[ 2 ] = 0;
			origin[ 3 ] = 0;

			region.SetSize( size );
			region.SetIndex( index );

			m_Output = OutputImageType::New();

			// set
			m_Output->SetRegions( region );
			m_Output->SetSpacing( spacing );
			m_Output->SetOrigin( origin );
			m_Output->Allocate();
			m_Output->FillBuffer( 0 );
		}


	}; // end class DkiSim

} // end namespace dki


/**
 * Constrained dkifit.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Diffusion kurtosis data simulator using two crossing fiber populations." );

	dki::DkiSim::sim_parameters args;

	args.sep = " ";

	args.size3D = std::vector< int >( 3 );
	args.size3D.at( 0 ) = 32;
	args.size3D.at( 1 ) = 32;
	args.size3D.at( 2 ) = 32;

	args.simBvals = std::vector< double >( 2 );
	args.simBvals.at( 0 ) = 1000;
	args.simBvals.at( 1 ) = 2000;

	args.numberOfBzeros = 1;
	args.numberOfOutliers = 1;
	args.B0Intensity = 100;
	args.SNR = std::numeric_limits< double >::infinity();
	args.FA = 0.7;
	args.MD = 0.0004;
	args.el1 = 0;
	args.el2 = 0;
	args.az1 = 0;
	args.az2 = 90;
	args.fraction1 = 0.5;
	args.fraction2 = 0.5;

	p.AddArgument( args.bvecsFileName, "bvecs" ) ->AddAlias( "r" ) ->SetDescription( "Diffusion gradient vectors (FSL format)" ) ->SetRequired( true );
	p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output tensor filename" ) ->SetRequired( true );
	p.AddArgument( args.sep, "separation" ) ->AddAlias( "s" ) ->SetDescription( "Separation string in bvecs file (default: ' ')" );

	p.AddArgument( args.size3D, "size" ) ->SetDescription( "Output 3D volume size (default: 32 32 32)" )->SetMinMax( 3, 3 );
	p.AddArgument( args.simBvals, "bvals" ) ->AddAlias( "b" ) ->SetDescription( "B-values to simulate (default: 1000 2000)" );
	p.AddArgument( args.numberOfBzeros, "number-of-bzeros" ) ->AddAlias( "nb" ) ->SetDescription( "Number of B-zero images (default: 1)" );
	p.AddArgument( args.numberOfOutliers, "number-of-outliers" ) ->AddAlias( "no" ) ->SetDescription( "Number of simulated outliers per voxel (default: 1)" );

	p.AddArgument( args.B0Intensity, "bzero-intensity" ) ->AddAlias( "bi" ) ->SetDescription( "B-zero intensity (default: 100)" );
	p.AddArgument( args.SNR, "signal-to-noise" ) ->AddAlias( "snr" ) ->SetDescription( "Signal to noise ratio (default: Inf)" );
	p.AddArgument( args.FA, "fractional-anisotropy" ) ->AddAlias( "fa" ) ->SetDescription( "Fractional anisotropy (default: 0.7)" );
	p.AddArgument( args.MD, "mean-diffusivity" ) ->AddAlias( "md" ) ->SetDescription( "Mean diffusivity (default: 0.0004)" );
	p.AddArgument( args.el1, "el1" ) ->AddAlias( "el1" ) ->SetDescription( "el1 (default: 0)" );
	p.AddArgument( args.el2, "el2" ) ->AddAlias( "el2" ) ->SetDescription( "el2 (default: 0)" );
	p.AddArgument( args.az1, "az1" ) ->AddAlias( "az1" ) ->SetDescription( "az1 (default: 0)" );
	p.AddArgument( args.az2, "az2" ) ->AddAlias( "az2" ) ->SetDescription( "az2 (default: 90)" );
	p.AddArgument( args.fraction1, "fiber-fraction1" ) ->AddAlias( "f1" ) ->SetDescription( "First fiber fraction (default: 0.5)" );
	p.AddArgument( args.fraction2, "fiber-fraction2" ) ->AddAlias( "f2" ) ->SetDescription( "Second fiber fraction (default: 0.5)" );


	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	dki::DkiSim sim( args );

	return EXIT_SUCCESS;
}
