#include "tkdCmdParser.h"

#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/accumulators/numeric/functional/complex.hpp>
#include <boost/accumulators/numeric/functional/valarray.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>

/**
 * Output max values and indices of given vector data in txt file.
 */
class AnnMax
{
public:

	typedef double ValueType;
	typedef std::vector< double > VectorType;
	typedef boost::tokenizer< boost::char_separator< char > > TokType;

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName, bool reportMin, bool reportMax, bool reportMinIndex, bool reportMaxIndex )
	{
		VectorType input = GetAverageVectorFromFile( inputFileName );

		ValueType min = GetMin( input );
		ValueType max = GetMax( input );
		int indexMin = GetIndex( input, min );
		int indexMax = GetIndex( input, max );

		if ( reportMin )
			std::cout << "min: " << min << std::endl;

		if ( reportMax )
			std::cout << "max: " << max << std::endl;

		if ( reportMinIndex )
			std::cout << "index min: " << indexMin << std::endl;

		if ( reportMaxIndex )
			std::cout << "index max: " << indexMax << std::endl;

		exit( EXIT_SUCCESS );
	}

protected:

	/**
	 * Return vector from file.
	 */
	VectorType GetVectorFromFile( const std::string& input )
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

		return output;
	}

	/**
	 * Return average vector from file.
	 */
	VectorType GetAverageVectorFromFile( const std::string& inputFileName )
	{

		std::ifstream in( inputFileName.c_str() );

		if ( !in )
		{
			std::cerr << "*** ERROR ***: Unable to open file \"" << inputFileName << "\" for input.\n";
			exit( EXIT_FAILURE );
		}

		std::string line;

		VectorType result;

		try
		{
			while ( getline( in, line ) )
			{
				ValueType value = 0;
				ValueType total = 0;
				TokType tok( line, boost::char_separator< char >( "\t" ) );

				for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
				{
					value += boost::lexical_cast< ValueType >( *id );
					total++;
				}

				// average...
				if ( total != 0 )
					value /= total;

				result.push_back( value );
			}
		}
		catch ( std::exception& e )
		{
			std::cerr << "*** ERROR ***: could not parse input. "
				"Are you sure columns are separated with the '\\t' char?" << std::endl;

			exit( EXIT_FAILURE );
		}

		in.close();

		return result;
	}

	/**
	 * Return max.
	 */
	ValueType GetMax( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::max > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return max( acc );
	}

	/**
	 * Return min.
	 */
	ValueType GetMin( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::min > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return min( acc );
	}

	/**
	 * Return index of given value (if not found return -1).
	 */
	int GetIndex( const VectorType& V, ValueType value )
	{
		int result = -1;
		for( unsigned int i = 0; i < V.size(); i++ )
		{
			if ( V[i] == value )
				return i;
		}
		return result;
	}

};

/**
 * Minimal connected components filter.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "ann max", "Output index and value of maximum value" );

	std::string inputFileName;
	bool reportMin = false;
	bool reportMax = false;
	bool reportMinIndex = false;
	bool reportMaxIndex = false;

	p.AddArgument( inputFileName, "input" )->AddAlias( "i" )->SetDescription( "Input file name" )->SetRequired( true );

	p.AddArgument( reportMin, "min" )->AddAlias( "l" )->SetDescription( "Print minimal value" )->SetRequired( false );
	p.AddArgument( reportMax, "max" )->AddAlias( "u" )->SetDescription( "Print maximal value" )->SetRequired( false );
	p.AddArgument( reportMinIndex, "min-index" )->AddAlias( "ml" )->SetDescription( "Print minimal value index (zero-based)" )->SetRequired( false );
	p.AddArgument( reportMaxIndex, "max-index" )->AddAlias( "mu" )->SetDescription( "Print maximal value index (zero-based)" )->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AnnMax annMax;
	annMax.Run( inputFileName, reportMin, reportMax, reportMinIndex, reportMaxIndex );

	return EXIT_SUCCESS;
}

