#include "tkdCmdParser.h"

#include "annCommon.h"

/**
 * Output positive and negative border distance.
 */
class AnnBorderDistance
{
public:

	typedef float PixelType;

	/**
	 * Run.
	 */
	void Run( const std::string& refFileName,
				const std::string& segFileName,
				const unsigned int minVoxels,
				const std::string& distanceMapFileName,
				const std::string& distanceMapInverseFileName,
				bool hausdorff )
	{

		PixelType negative;
		PixelType positive;
		ann::Ann< PixelType >::BorderDistance( positive, negative, refFileName, segFileName,
				minVoxels, distanceMapFileName, distanceMapInverseFileName );

		std::cout << "False positive border distance: " << positive << std::endl;
		std::cout << "False negative border distance: " << negative << std::endl;

		if ( hausdorff )
		{
			PixelType negative;
			PixelType positive;

			ann::Ann< PixelType >::HausdorffDistance( negative, positive, refFileName, segFileName,
					minVoxels );

			std::cout << "False negative Hausdorff distance: " << negative << std::endl;
			std::cout << "False positive Hausdorff distance: " << positive << std::endl;
		}
	}
};

/**
 * Minimal connected components filter.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann border distance", "Output average positive and negative border distance for minimal connected components" );

	std::string refFileName;
	std::string segFileName;

	int minVoxels = 0;
	std::string distanceMapFileName;
	std::string distanceMapInverseFileName;
	bool hausdorff;

	p.AddArgument( refFileName, "ref" )->AddAlias( "r" )->SetDescription( "Reference image file" )->SetRequired( true );
	p.AddArgument( segFileName, "seg" )->AddAlias( "s" )->SetDescription( "Segmented image file" )->SetRequired( true );
	p.AddArgument( minVoxels, "min-voxels" )->AddAlias( "min" )->SetDescription( "Minimal voxels in isolated components (default: 0)" )->SetRequired( false );

	p.AddArgument( distanceMapFileName, "distance-map" )->AddAlias( "d" )->SetDescription( "Positive distance map image file (DEBUG; default: empty)" )->SetRequired( false );
	p.AddArgument( distanceMapInverseFileName, "distance-map-inverse" )->AddAlias( "di" )->SetDescription( "Negative distance map image file (DEBUG; default: empty)" )->SetRequired( false );
	p.AddArgument( hausdorff, "hausdorff" )->AddAlias( "h" )->SetDescription( "Hausdorff distance (default: false)" )->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnBorderDistance borderDistance;
	borderDistance.Run( refFileName, segFileName, minVoxels, distanceMapFileName, distanceMapInverseFileName, hausdorff );

	return EXIT_SUCCESS;
}


