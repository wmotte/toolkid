#include "tkdCmdParser.h"
#include "itkImage.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <boost/random.hpp>
//#include <boost/test/unit_test.hpp>
//#include <boost/test/floating_point_comparison.hpp>

#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/accumulators/numeric/functional/complex.hpp>
#include <boost/accumulators/numeric/functional/valarray.hpp>

#include <boost/accumulators/accumulators.hpp>

#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/p_square_cumulative_distribution.hpp>

#include "graphCommon.h"

using namespace boost;
using namespace boost::accumulators;


/**
 * @author: W.M. Otte
 * @date: 19-11-2009
 * @aim: Determine histogram from given input numbers.
 * To determine image histograms instead of list histograms, use imagehistogram.
 */
class Histogram
{

public:

	typedef double PixelType;
	typedef vnl_vector< PixelType > VectorType;
	typedef vnl_matrix< PixelType > MatrixType;

	/**
	 * Run.
	 */
	void run( const std::string& input, const std::string& output,
					bool cumulative, bool reverse, unsigned int bins )
	{
		// [ 1 ] Read time-series.
		// --------------------------------------
		VectorType inputs;
		graph::Graph< PixelType >::GetVectorFromFile( input, inputs );

		// [ 2 ] get histogram.
		// --------------------------------------
		VectorType positions( bins );
		VectorType values( bins );

		getHistogram( positions, values, inputs, bins, cumulative, reverse );

		// [ 3 ] write to disk.
		// --------------------------------------

		graph::Graph< PixelType >::WriteVectorToFile( output + "_density.txt", values );

		graph::Graph< PixelType >::WriteVectorToFile( output + "_positions.txt", positions );

		exit( EXIT_SUCCESS );
	}

protected:

	/**
	 * Convert input vector into histogram.
	 *
	 * The histogram density estimator returns a histogram of the sample distribution.
	 * The positions and sizes of the bins are determined
	 * using a specifiable number of cached samples (cache_size).
	 *
	 * The range between the minimum and the maximum of the cached samples
	 * is subdivided into a specifiable number of bins (num_bins) of same size.
	 * Additionally, an under- and an overflow bin is added to
	 * capture future under- and overflow samples.
	 * Once the bins are determined, the cached samples and all
	 * subsequent samples are added to the correct bins.
	 * At the end, a range of std::pair is return,
	 * where each pair contains the position of the bin (lower bound)
	 * and the samples count (normalized with the total number of samples).
	 *
	 * density_cache_size: Number of first samples used to determine min and max.
	 * density_num_bins: Number of bins (two additional bins collect under- and overflow samples).
	 */
	void getHistogram( VectorType& positions, VectorType& values,
							const VectorType& input, unsigned int bins,
								bool cumulative, bool reverse )
	{
		typedef accumulator_set< PixelType, stats< tag::density > > acc_density_t;
		typedef accumulator_set< PixelType, stats< tag::p_square_cumulative_distribution > > acc_cumul_t;
		typedef iterator_range< std::vector< std::pair< PixelType, PixelType > >::iterator > iterator_type;

		if ( ! cumulative )
		{
			// define accumulator with density feature.
			acc_density_t acc( tag::density::cache_size = input.size(), tag::density::num_bins = bins );

			// fill accumulator
			for ( unsigned int i = 0; i < input.size(); i++ )
			{
				acc( input[i] );
			}

			// get densities
			iterator_type densities = density( acc );

			positions = VectorType( densities.size() );
			values    = VectorType( densities.size() );

			// insert densities in output vector.
			for( unsigned int i = 0; i < densities.size(); i++ )
			{
				positions[i] = densities[i].first;
				values[i]    = densities[i].second;
			}

		}
		else // cumulative
		{
			acc_cumul_t acc( p_square_cumulative_distribution_num_cells = bins );

			// fill accumulator
			for ( unsigned int i = 0; i < input.size(); i++ )
			{
				acc( input[i] );
			}
			iterator_type histogram = p_square_cumulative_distribution( acc );

			// insert histogram in output vector.
			for( unsigned int i = 0; i < histogram.size(); i++ )
			{
				positions[i] = i;
				values[i]    = histogram[i].first;
			}

		}

		if ( reverse )
		{
			values.flip();
		}
	}
};

/**
 * Main hierarchy.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "histogram", "Construct (cumulative) histogram from given list of numbers." );

	std::string input;
	std::string output;

	bool cumulative = false;
	bool reverse = false;
	int bins = 128;

	p.AddArgument( input, "input" )
		->AddAlias( "i" )
		->SetInput( "filename" )
		->SetDescription( "List input (column format)" )
		->SetRequired( true )
		-> SetMinMax( 1, 1 );

	p.AddArgument( output, "output" )
		->AddAlias( "o" )
		->SetInput( "filename" )
		->SetDescription( "Rootname histogram output (column format)" )
		->SetRequired( true )
		-> SetMinMax( 1, 1 );

	p.AddArgument( cumulative, "cumulative" )
		->AddAlias( "c" )
		->SetInput( "bool" )
		->SetDescription( "Cumulative histogram [ default: false ]" )
		->SetRequired( false );

	p.AddArgument( reverse, "reverse" )
		->AddAlias( "r" )
		->SetInput( "bool" )
		->SetDescription( "Reverse histogram [ default: false ]" )
		->SetRequired( false );

	p.AddArgument( bins, "bins" )
		->AddAlias( "b" )
		->SetInput( "uint" )
		->SetDescription( "Number of histogram bins [ default: 128 ]" )
		->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Histogram histogram;

	histogram.run( input, output, cumulative, reverse, ( unsigned ) bins );

	return EXIT_SUCCESS;
}
