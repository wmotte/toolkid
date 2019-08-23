#include "powerlawCommon.h"
//#include <boost/program_options.hpp>
#include "tkdCmdParser.h"

/**
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC Utrecht, NL.
 * @date: 19-11-2009
 *
 * Estimate powerlaw scaling parameter from input distribution.
 *
 * ***************************************************************************
 * Method: "Power-law distributions in empirical data", Clauset et al, 2009
 * http://www.santafe.edu/~aaronc/powerlaws/
 * ***************************************************************************
 */
class PowerLawFit
{

public:

	typedef double ValueType;
	typedef std::vector< ValueType > VectorType;

	/**
	 * Power law fit.
	 */
	void run( const std::string& inputFileName, bool nosmall, bool finite,
					double startXmin, double incrementXmin, double endXmin,
						bool bootstrap, unsigned int bootstrapIterations, bool verbose )
	{
		// [ 1 ] read input from text file ...
		VectorType values = getInput( inputFileName );

		// [ 2 ] bootstrap or single fit ...
		VectorType results;

		if ( bootstrap )
		{
			graph::Powerlaw< ValueType >::BootstrapFit( values, results, nosmall, finite, startXmin, incrementXmin, endXmin, bootstrapIterations, verbose );

			if ( ! results.empty() )
			{
				std::cout << "Alpha," << results.at( 0 ) <<  std::endl;
				std::cout << "Xmin," << results.at( 1 ) << std::endl;
				std::cout << "Log-likelihood," << results.at( 2 ) << std::endl;
				std::cout << "Alpha_sd," << results.at( 3 ) <<  std::endl;
				std::cout << "Xmin_sd," << results.at( 4 ) << std::endl;
				std::cout << "Log-likelihood_sd," << results.at( 5 ) << std::endl;

			}
			else
			{
				std::cerr << "*** ERROR ***: maximum likelihood "
					"bootstrap estimation failed! -> check input ..." << std::endl;
			}

		}
		else
		{
			graph::Powerlaw< ValueType >::SingleFit( values, results, nosmall, finite,
					startXmin, incrementXmin, endXmin );

			if ( ! results.empty() )
			{
				std::cout << "Alpha," << results.at( 0 ) << std::endl;
				std::cout << "Xmin," << results.at( 1 ) << std::endl;
				std::cout << "Log-likelihood," << results.at( 2 ) << std::endl;
			}
			else
			{
				std::cerr << "*** ERROR ***: maximum likelihood "
					"single estimation failed! -> check input ..." << std::endl;
			}
		}
	}

protected:

	/**
	 * Return input from given text file as vector.
	 */
	VectorType getInput( const std::string& input )
	{
		std::ifstream inFile;

		inFile.open( input.c_str() );
		if ( !inFile )
		{
			std::cout << "*** ERROR ***: Unable to open: " << input << "." << std::endl;
			exit( EXIT_FAILURE );
		}

		double x;
		VectorType output;

		while ( inFile >> x )
		{
			output.push_back( x );
		}
		inFile.close();

		/**
		 * Negative values will be converted to complex numbers in matlab,
		 * but not with the stl ...
		 *
		 * No support is given (yet) for complex number mle.
		 */
		if ( *( std::min_element( output.begin(), output.end() ) ) < 0 )
		{
			std::cerr << "*** ERROR ***: Negative input not supported!" << std::endl;
			exit (EXIT_FAILURE );
		}

		return output;
	}
};

// ************************************************************************************


/**
 * Option list.
 */
struct parameters
{
	std::string input;
	bool nosmall;
	bool finite;
	bool bootstrap;
	bool verbose;
	double startXmin;
	double incrementXmin;
	double endXmin;
	int bootstrapIterations;

};

/**
 * Fit powerlaw to list of numbers.
 */
int main(int argc, char* argv[])
{
	tkd::CmdParser p( argv[0], "Fits a power-law distributional model to data." );

		parameters list;

		list.nosmall = false;
		list.finite = false;
		list.bootstrap = false;
		list.verbose = false;
		list.startXmin = 1.5;
		list.incrementXmin = 0.01;
		list.endXmin = 3.5;
		list.bootstrapIterations = 1000;

		p.AddArgument( list.input, "input" )
				->AddAlias( "i" )
				->SetDescription( "Input file with distribution values in column format" )
				->SetRequired( true );

		p.AddArgument( list.nosmall, "nosmall" )
				->SetDescription( "Truncate the search over xmin values before the finite-size bias becomes significant (default: false)" );

		p.AddArgument( list.finite, "finite" )
				->SetDescription( "Use an experimental finite-size correction (default: false)" );

		p.AddArgument( list.bootstrap, "bootstrap" )
				->SetDescription( "Run non-parametric bootstrap instead of single estimation (default: false)" );

		p.AddArgument( list.verbose, "verbose" )
				->SetDescription( "Print boostrap status (default: false)" );

		p.AddArgument( list.startXmin, "start-x-min" )
				->SetDescription( "Start value for discrete xmin estimation (default: 1.5)" );

		p.AddArgument( list.incrementXmin, "increment-x-min" )
				->SetDescription( "Increment value for discrete xmin estimation (default: 0.01)" );

		p.AddArgument( list.endXmin, "end-x-min" )
				->SetDescription( "End value for discrete xmin estimation (default: 3.5)" );

		p.AddArgument( list.bootstrapIterations, "iterations" )
				->SetDescription( "Bootstrap iterations (default: 1000)" );


		if ( !p.Parse( argc, argv ) )
		{
			p.PrintUsage( std::cout );
			return EXIT_FAILURE;
		}

        // run application ...
        PowerLawFit powerlawFit;

    	powerlawFit.run( list.input, list.nosmall, list.finite,
							list.startXmin, list.incrementXmin, list.endXmin,
											list.bootstrap, list.bootstrapIterations, list.verbose );

    return EXIT_SUCCESS;
}
