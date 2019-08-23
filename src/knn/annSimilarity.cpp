#include "tkdCmdParser.h"

#include "annCommon.h"

/**
 * Calculate Similarity index and probabilistic similarity index (anbeek, neuroimage).
 * Added: specificity and sensitivity.
 */
class AnnSimilarity
{
public:

	typedef float PixelType;

	/**
	 * Run.
	 */
	void Run( const std::string& ref,
				const std::string& seg,
				bool si,
				bool probabilistic,
				bool specificity,
				bool sensitivity,
				const std::string& maskFileName )
	{
		if ( si )
		{
			std::cout << "SI," << ann::Ann< PixelType >::SI( ref, seg ) << std::endl;
		}

		// psi
		if ( probabilistic )
		{
			std::cout << "PSI," << ann::Ann< PixelType >::PSI( ref, seg ) << std::endl;
		}

		if ( sensitivity )
		{
			std::cout << "Sensitivity," << ann::Ann< PixelType >::Sensitivity( ref, seg ) << std::endl;
		}

		if ( specificity )
		{
			if ( !maskFileName.empty() )
			{
				std::cout << "Specificity," << ann::Ann< PixelType >::Specificity( ref, seg, maskFileName ) << std::endl;
			}
			else
			{
				std::cerr << "*** ERROR ***: no mask supplied!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}
	}
};

/**
 * Similarity index.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann similarity", "Calculate (probabilistic) similarity index" );

	std::string ref;
	std::string seg;

	std::string maskFileName;

	bool si = false;
	bool probabilisitic = false;
	bool specificity = false;
	bool sensitivity = false;

	p.AddArgument( ref, "ref" )->AddAlias( "r" )->SetDescription( "Reference input image" )->SetRequired( true );
	p.AddArgument( seg, "seg" )->AddAlias( "s" )->SetDescription( "(Prob) Segmentation input image" )->SetRequired( true );

	p.AddArgument( probabilisitic, "probabilistic" )->AddAlias( "psi" )->SetDescription( "Run probabilistic SI [default: false]" )->SetRequired( false );
	p.AddArgument( si, "si" )->AddAlias( "si" )->SetDescription( "Run SI [default: false]" )->SetRequired( false );
	p.AddArgument( maskFileName, "mask" )->AddAlias( "m" )->SetDescription( "Mask input image" )->SetRequired( false );
	p.AddArgument( specificity, "specificity" )->AddAlias( "sp" )->SetDescription( "Output specificity [default: false]" )->SetRequired( false );
	p.AddArgument( sensitivity, "sensitivity" )->AddAlias( "se" )->SetDescription( "Output sensitivity [default: false]" )->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnSimilarity annSimilarity;
	annSimilarity.Run( ref, seg, si, probabilisitic, specificity, sensitivity, maskFileName );

	return EXIT_SUCCESS;
}


