#include "tkdCmdParser.h"

#include "annCommon.h"

/**
 * Output Stats segmentation data.
 */
class AnnStats
{
public:

	typedef float PixelType;

	/**
	 * Run.
	 */
	void Run( const std::string& ref,
				const std::string& seg,
				const std::string& mask,
				const std::string& outputSen,
				const std::string& outputSpec,
				const std::string& outputAuc,
				const std::string& outputSI,
				const std::string& outputPositiveBorderDistance,
				const std::string& outputNegativeBorderDistance,
				const std::string& outputPositiveHaussdorff,
				const std::string& outputNegativeHaussdorff,
				const unsigned int minVoxels )
	{
		if ( !outputPositiveHaussdorff.empty() || !outputNegativeHaussdorff.empty() )
		{
			ann::Ann< PixelType >::VectorType negative;
			ann::Ann< PixelType >::VectorType positive;
			ann::Ann< PixelType >::HausdorffDistanceRange( negative, positive, ref, seg, minVoxels );

			if ( !outputNegativeHaussdorff.empty() )
				ann::Ann< PixelType >::WriteData( negative, outputNegativeHaussdorff );

			if ( !outputPositiveHaussdorff.empty() )
					ann::Ann< PixelType >::WriteData( positive, outputPositiveHaussdorff );
		}

		if( !outputNegativeBorderDistance.empty() || !outputPositiveBorderDistance.empty() )
		{
			ann::Ann< PixelType >::VectorType negative;
			ann::Ann< PixelType >::VectorType positive;
			ann::Ann< PixelType >::BorderDistanceRange( negative, positive, ref, seg, minVoxels );

			if ( !outputNegativeBorderDistance.empty() )
				ann::Ann< PixelType >::WriteData( negative, outputNegativeBorderDistance );

			if ( !outputPositiveBorderDistance.empty() )
					ann::Ann< PixelType >::WriteData( positive, outputPositiveBorderDistance );
		}

		if ( !outputAuc.empty() || !outputSen.empty() || !outputSpec.empty() )
		{
			ann::Ann< PixelType >::VectorType sensitivity;
			ann::Ann< PixelType >::VectorType specificity; // ~ 1 - specificity!
			ann::Ann< PixelType >::VectorType roc;

			PixelType auc = ann::Ann< PixelType >::ROC( sensitivity, specificity, roc, ref, seg, mask, minVoxels );
			std::cout << "AUC," << auc << std::endl;

			if ( !outputSen.empty() )
				ann::Ann< PixelType >::WriteData( sensitivity, outputSen );

			if ( !outputSpec.empty() )
				ann::Ann< PixelType >::WriteData( specificity, outputSpec );

			if ( !outputAuc.empty() )
				ann::Ann< PixelType >::WriteData( roc, outputAuc );
		}

		if ( !outputSI.empty() )
		{
			ann::Ann< PixelType >::VectorType si;
			ann::Ann< PixelType >::RangeSI( si, ref, seg, minVoxels );
			ann::Ann< PixelType >::WriteData( si, outputSI );
		}
	}
};

/**
 * Stats.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann stats", "Calculate segmentation statistics." );

	std::string ref;
	std::string seg;
	std::string outputSen;
	std::string outputSpec;
	std::string outputAuc;
	std::string outputSI;
	std::string maskFileName;
	std::string outputPositiveBorderDistance;
	std::string outputNegativeBorderDistance;
	std::string outputNegativeHausdorff;
	std::string outputPositiveHausdorff;

	int minVoxels = 0;

	p.AddArgument( ref, "ref" )->AddAlias( "r" )->SetDescription( "Reference input image" )->SetRequired( true );
	p.AddArgument( seg, "seg" )->AddAlias( "s" )->SetDescription( "(Prob) Segmentation input image" )->SetRequired( true );
	p.AddArgument( maskFileName, "mask" )->AddAlias( "m" )->SetDescription( "Mask input image" )->SetRequired( true );
	p.AddArgument( outputSen, "output-sensitivity" )->AddAlias( "ose" )->SetDescription( "Output sensitivity file" )->SetRequired( false );
	p.AddArgument( outputSpec, "output-specificity" )->AddAlias( "osp" )->SetDescription( "Output specificity file (~ 1 - specificity)" )->SetRequired( false );
	p.AddArgument( outputAuc, "output-auc" )->AddAlias( "oa" )->SetDescription( "Output cumulative AUC file" )->SetRequired( false );
	p.AddArgument( outputSI, "output-si" )->AddAlias( "si" )->SetDescription( "Output similarity index file" )->SetRequired( false );
	p.AddArgument( outputPositiveBorderDistance, "output-positive-border-distance" )->AddAlias( "opd" )->SetDescription( "Output positive border distance file" )->SetRequired( false );
	p.AddArgument( outputNegativeBorderDistance, "output-negative-border-distance" )->AddAlias( "ond" )->SetDescription( "Output negative border distance file" )->SetRequired( false );
	p.AddArgument( outputPositiveHausdorff, "output-positive-hausdorff-distance" )->AddAlias( "oph" )->SetDescription( "Output positive Hausdorff distance file" )->SetRequired( false );
	p.AddArgument( outputNegativeHausdorff, "output-negative-hausdorff-distance" )->AddAlias( "onh" )->SetDescription( "Output negative Hausdorff distance file" )->SetRequired( false );
	p.AddArgument( minVoxels, "min-voxels" )->AddAlias( "min" )->SetDescription( "Minimal voxels in isolated components (default: disabled, enable with non-zero input)" )->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnStats annStats;
	annStats.Run( ref, seg, maskFileName, outputSen, outputSpec, outputAuc, outputSI, outputPositiveBorderDistance, outputNegativeBorderDistance,
			outputPositiveHausdorff, outputNegativeHausdorff, (unsigned) minVoxels );

	return EXIT_SUCCESS;
}


