#include "tkdCmdParser.h"

#include "graphCommon.h"
#include <limits>
#include <stdlib.h> // RAND_MAX def.
#include <stdio.h>
#include <time.h>

/**
 * Option list.
 */
struct parameters
{
	std::string inputFileName;
	std::string outputFileName;
	float lambda;
	int number_iterations;
	int number_of_exemplars;
	float preference;
	bool removeDegeneracies;
};

/**
 * Affinity propagation clustering.
 *
 * Affinity propagation( Science vol 315, 2007, Frey and Dueck).
 */
class AffinityPropagation
{
public:

	typedef double PixelType;

	/**
	 * Run app.
	 */
	void Run( parameters& list )
	{
		if ( graph::Graph< PixelType >::GetImageDimensions( list.inputFileName ) == 2 )
		{
			graph::Graph< PixelType >::ImagePointerType image;
			graph::Graph< PixelType >::ReadMatrix( image, list.inputFileName );

			PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
			graph::Graph< PixelType >::ImageType::RegionType region = image->GetLargestPossibleRegion();
			graph::Graph< PixelType >::ImageType::SizeType size = region.GetSize();

			int rows = size[0];
			int cols = size[1];
			graph::Graph< PixelType >::DataMatrixType G( rows, cols, buffer );

			// similarity matrix S
			graph::Graph< PixelType >::MatrixType S( rows, cols );

			for ( int i = 0; i < rows; i++ )
				for ( int j = 0; j < cols; j++ )
					S( i, j ) = G( i, j );

			if( list.preference > std::numeric_limits< float >::min() )
			{
				S.fill_diagonal( list.preference );
			}
			else
			{
				// get median of input matrix.
				graph::Graph< PixelType >::VectorType medians( S.rows(), 0 );
				for( unsigned int i = 0; i < S.rows(); i++ )
					medians( i ) = graph::Graph< PixelType >::GetMedian( S.get_row( i ) );

				// set all preferences to minimum median.
				PixelType min = graph::Graph< PixelType >::GetMin( medians );
				S.fill_diagonal( min );
				std::cout << "Minimal median: " << min << std::endl;
			}

			// cluster
			graph::Graph< PixelType >::VectorType clusters;
			AP( clusters, S, list );

			// Write
			graph::Graph< PixelType >::WriteVectorToFile( list.outputFileName, clusters );
		} else
		{
			std::cerr << "Number of input (matrix) dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	/**
	 * Affinity propagation of given similarity matrix.
	 */
	void AP( graph::Graph< PixelType >::VectorType& clusters, graph::Graph< PixelType >::MatrixType& Similarities, parameters& list )
	{

		int number_of_nodes = Similarities.rows();

		// get min and max.
		PixelType min = std::numeric_limits< PixelType >::max();
		PixelType max = std::numeric_limits< PixelType >::min();
		for ( unsigned int i = 0; i < Similarities.rows(); i++ )
		{
			for ( unsigned j = 0; j < Similarities.cols(); j++ )
			{
				if ( Similarities( i, j ) < min )
					min = Similarities( i, j );
				if ( Similarities( i, j ) > max )
					max = Similarities( i, j );
			}
		}

		if( list.removeDegeneracies )
		{
			std::cout << "Adding noise..." << std::endl;
			srand( time( NULL ) );
			PixelType diff = max - min;

			// add noise to similarity matrix.
			for ( unsigned int i = 0; i < Similarities.rows(); i++ )
			{
				for ( unsigned int j = 0; j < Similarities.cols(); j++ )
				{
					if( ( static_cast< PixelType >( rand() ) / RAND_MAX ) < 0.5 )
						Similarities( i, j ) += 1.0E-12 * ( - static_cast< PixelType >( rand() ) / RAND_MAX ) * diff;
					else
						Similarities( i, j ) += 1.0E-12 * ( static_cast< PixelType >( rand() ) / RAND_MAX ) * diff;
				}
			}
		}

		// initialize availability and responsibility graphs.
		graph::Graph< PixelType >::MatrixType Availabilities( number_of_nodes, number_of_nodes, 0 );
		graph::Graph< PixelType >::MatrixType Responsibilities( number_of_nodes, number_of_nodes, 0 );

		// iterate.
		for ( int iteration = 0; iteration < list.number_iterations; iteration++ )
		{
			std::cerr << "Iteration: " << iteration << std::endl;

			graph::Graph< PixelType >::MatrixType ROld( Responsibilities );
			graph::Graph< PixelType >::MatrixType AS = Similarities + Availabilities;

			// max AS...
			std::vector< PixelType > max_values( number_of_nodes, 0 );
			std::vector< int > index_of_max_values( number_of_nodes, 0 );

			for ( int i = 0; i < number_of_nodes; i++ )
			{
				PixelType maximum = std::numeric_limits< PixelType >::min();
				int position = -1;
				for ( int j = 0; j < number_of_nodes; j++ )
				{
					if ( AS( i, j ) <= maximum )
						continue;

					maximum = AS( i, j );
					position = j;
				}

				max_values.at( i ) = maximum;
				index_of_max_values.at( i ) = position;
			}

			for ( int i = 0; i < number_of_nodes; i++ )
			{
				AS( i, index_of_max_values.at( i ) ) = -std::numeric_limits< PixelType >::infinity();
			}

			// max AS second time...
			std::vector< PixelType > max_values2( number_of_nodes, 0 );
			std::vector< int > index_of_max_values2( number_of_nodes, 0 );

			for ( int i = 0; i < number_of_nodes; i++ )
			{
				PixelType maximum = std::numeric_limits< PixelType >::min();
				int position = -1;
				for ( int j = 0; j < number_of_nodes; j++ )
				{
					if ( AS( i, j ) <= maximum )
						continue;

					maximum = AS( i, j );
					position = j;
				}

				max_values2.at( i ) = maximum;
				index_of_max_values2.at( i ) = position;
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					Responsibilities( i, j ) = Similarities( i, j ) - max_values.at( i );
				}
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				Responsibilities( i, index_of_max_values.at( i ) ) = Similarities( i, index_of_max_values.at( i ) ) - max_values2.at( i );
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					Responsibilities( i, j ) = ( 1.0 - list.lambda ) * Responsibilities( i, j ) + list.lambda * ROld( i, j );
				}
			}

			graph::Graph< PixelType >::MatrixType AOld( Availabilities );
			graph::Graph< PixelType >::MatrixType R_positive_values( number_of_nodes, number_of_nodes, 0 );

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					if ( Responsibilities( i, j ) >= 0.0 )
					{
						R_positive_values( i, j ) = Responsibilities( i, j );
					}
				}
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				R_positive_values( i, i ) = Responsibilities( i, i );
			}

			std::vector< PixelType > sum_matrix( number_of_nodes, 0 );

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					sum_matrix.at( i ) += R_positive_values( j, i );
				}
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					Availabilities( j, i ) = sum_matrix.at( i ) - R_positive_values( j, i );
				}
			}

			std::vector< PixelType > diag_matrix( number_of_nodes, 0 );
			for ( int i = 0; i < number_of_nodes; ++i )
			{
				diag_matrix.at( i ) = Availabilities( i, i );
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					if ( Availabilities( i, j ) >= 0.0 )
					{
						Availabilities( i, j ) = 0.0;
					}
				}
			}

			for ( int i = 0; i < number_of_nodes; ++i )
			{
				Availabilities( i, i ) = diag_matrix.at( i );
			}
			for ( int i = 0; i < number_of_nodes; ++i )
			{
				for ( int j = 0; j < number_of_nodes; ++j )
				{
					Availabilities( i, j ) = ( 1.0 - list.lambda ) * Availabilities( i, j ) + list.lambda * AOld( i, j );
				}
			}
		}

		// Pseudomarginals
		graph::Graph< PixelType >::MatrixType Evidence( number_of_nodes, number_of_nodes, 0 );
		Evidence = Responsibilities + Availabilities;

		// TODO get diagonal. Values above zero are exemplars.
		for( int i = 0; i < number_of_nodes; i++ )
			std::cout << Evidence( i, i ) << std::endl;

		std::cerr << "End evidence..." << std::endl;

		//TODO
		//std::cout << "Writing evidence matrix..." << std::endl;
		//graph::Graph< PixelType >::WriteMatrixToFile( Evidence, "evidence.nii.gz" );


		std::vector< int > pos_diag_idx;
		for ( int i = 0; i < number_of_nodes; ++i )
		{
			if ( Evidence( i, i ) > 0.0 )
			{
				pos_diag_idx.push_back( i );
			}
		}

		list.number_of_exemplars = pos_diag_idx.size();
		if ( list.number_of_exemplars == 0 )
		{
			list.number_of_exemplars = 1;
			pos_diag_idx.push_back( 0 );
		}

		graph::Graph< PixelType >::MatrixType exemplar_cols( number_of_nodes, list.number_of_exemplars );
		for ( int i = 0; i < number_of_nodes; ++i )
		{
			for ( int j = 0; j < list.number_of_exemplars; ++j )
			{
				exemplar_cols( i, j ) = Similarities( i, pos_diag_idx.at( j ) );
			}
		}

		std::vector< PixelType > tmp( number_of_nodes );
		std::vector< int > c( number_of_nodes );

		for ( int i = 0; i < number_of_nodes; ++i )
		{
			PixelType maximum = std::numeric_limits< PixelType >::min();
			int position = -1;
			for ( int j = 0; j < list.number_of_exemplars; ++j )
			{
				if ( exemplar_cols( i, j ) <= maximum )
					continue;
				maximum = exemplar_cols( i, j );
				position = j;
			}

			tmp.at( i ) = maximum;
			c.at( i ) = position;
		}

		for ( unsigned int i = 0; i < pos_diag_idx.size(); ++i )
		{
			c.at( pos_diag_idx.at( i ) ) = i;
		}

		std::vector< int > idx( number_of_nodes );
		for ( int i = 0; i < number_of_nodes; ++i )
		{
			idx.at( i ) = pos_diag_idx.at( c.at( i ) );
		}


		//for( unsigned int i = 0; i < c.size(); i++ )
		//{
		//	std::cout << c.at( i ) << std::endl;
		//}

		//clusters = graph::Graph< PixelType >::VectorType( number_of_nodes );
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Affinity propagation clustering of input similarity matrix" );

	parameters list;
	list.lambda = 5.;
	list.number_iterations = 100;
	list.number_of_exemplars = 0;
	list.preference = std::numeric_limits< float >::min();
	list.removeDegeneracies = false;

	p.AddArgument( list.inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription(
			"Input image: 2D adjacency matrix" ) ->SetRequired( true );

	p.AddArgument( list.outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription(
			"Output file: cluster txt file" ) ->SetRequired( true );

	p.AddArgument( list.lambda, "lambda" ) ->AddAlias( "l" ) ->SetInput( "float" ) ->SetDescription(
			"Damping factor (dafault: 0.5)" ) ->SetRequired( false );

	p.AddArgument( list.preference, "preference" ) ->AddAlias( "p" ) ->SetInput( "float" ) ->SetDescription(
			"Preference (default: minimum similarity of the median similarity)" ) ->SetRequired( false );

	p.AddArgument( list.number_iterations, "iterations" ) ->AddAlias( "it" ) ->SetInput( "int" ) ->SetDescription(
			"Iterations (default: 100)" ) ->SetRequired( false );

	p.AddArgument( list.removeDegeneracies, "remove-degeneracies" ) ->AddAlias( "d" ) ->SetInput( "bool" ) ->SetDescription(
			"Remove degeneracies by adding small noise to the data (default: false)" ) ->SetRequired( false );



	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	AffinityPropagation ap;
	ap.Run( list );

	return EXIT_SUCCESS;
}
