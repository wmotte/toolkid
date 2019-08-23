#include "tkdCmdParser.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/algo/vnl_svd.h"

#include <stdlib.h>
#include <math.h>
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/utility.hpp>
#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/accumulators/numeric/functional/complex.hpp>
#include <boost/accumulators/numeric/functional/valarray.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <boost/math/distributions/normal.hpp>

/* Indexer problem... */

/**
 * Option list.
 */
struct parameters
{
	std::string inputFileName;
	std::string link;
	std::string separation;
	int iterations;
	double percentage;
	bool verbose;
	bool printPValues;
	bool printZValues;
};

/**
 * Generalized Linear Model (GLM) is a flexible generalization of ordinary least squares regression.
 *
 * http://en.wikipedia.org/wiki/Generalized_linear_model
 *
 * wim@invivonmr.uu.nl, 02-03-2010.
 *
 */
namespace generalized_linear_model
{
	typedef double PixelType;
	typedef vnl_matrix< PixelType > MatrixType;
	typedef vnl_vector< PixelType > VectorType;
	typedef boost::tokenizer< boost::char_separator< char > > TokType;
	typedef boost::mt19937 random_number_type;
	typedef boost::uniform_int< long > int_distribution_type;
	typedef boost::variate_generator< random_number_type&, int_distribution_type > int_generator_type;

	class GLM
	{
	public:

		/**
		 * Run glm.
		 */
		void Run( parameters& list )
		{
			// read covariates and response-function (latter as last column).
			MatrixType covariates;
			VectorType binaryResponse;

			if( list.verbose )
				std::cout << "Reading input data..." << std::endl;

			ReadData( list.separation, list.inputFileName, covariates, binaryResponse, list.percentage );
			covariates = AddIntercept( covariates );

			if( list.verbose )
				std::cout << "Finished reading input data..." << std::endl;

			if ( list.link == "logistic" )
			{
				// Calculate coefficients...
				VectorType coefficients = IRLS( binaryResponse, covariates, list.iterations, list.verbose );
				std::cout << "Coefficients: " << coefficients << std::endl;

				if( list.printZValues )
					std::cout << "Z-values: " << ZStatistic( coefficients, binaryResponse, covariates ) << std::endl;

				if( list.printPValues )
					std::cout << "P-values: " << PValue( ZStatistic( coefficients, binaryResponse, covariates ) ) << std::endl;
			}
			else
			{
				std::cerr << "*** ERROR ***: Link function not supported!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

	protected:

		/**
		 * The P values for individual coefficients, given the z values.
		 */
		VectorType PValue( const VectorType& zValues )
		{
			VectorType pValues( zValues.size(), 0.0 );

			for ( unsigned int i = 0; i < zValues.size(); i++ )
			{
				// two-sided
				boost::math::normal s;
				PixelType alpha = boost::math::pdf( s, zValues( i ) ) * 2.;
				pValues( i ) = alpha;
			}
			return pValues;
		}

		/**
		 * The z statistics for individual coefficients.
		 */
		VectorType ZStatistic( const VectorType& coefficients, const VectorType& response,
				const MatrixType& covariates )
		{
			VectorType means = Means( coefficients, covariates );
			MatrixType weights = Weights( means, coefficients, covariates );

			MatrixType xwx = covariates * weights * covariates.transpose();

			vnl_svd< double > svd( xwx );

			MatrixType variance = svd.inverse();

			VectorType coefficientSE( variance.rows(), 0.0 );

			for( unsigned int i = 0; i < variance.rows(); i++ )
				coefficientSE( i ) = std::sqrt( variance( i, i ) );

			MatrixType correlation( variance.rows(), variance.rows(), 0.0 );

			for( unsigned int i = 0; i < variance.rows(); i++)
			{
				for ( unsigned int j = i; j < variance.rows(); j++ )
				{
					correlation( i, j ) = variance( i, j ) /
											std::sqrt(variance( i, i ) * variance( j, j ) );

					correlation( j, i ) = correlation( i, j );
				}
			}

			VectorType zValues( coefficients.size(), 0.0 );

			for( unsigned int i = 0; i < zValues.size(); i++ )
				zValues( i ) = coefficients( i ) / coefficientSE( i );

			return zValues;
		}

		/**
		 * Return number of training dimensions.
		 */
		unsigned int GetNumberOfTrainingDims( const std::string& trainingData )
		{
			std::string sep = "\t";
			std::ifstream in( trainingData.c_str() );

			if ( in.fail() )
			{
				std::cerr << "*** ERROR ***: could not read labels from: " << trainingData << std::endl;
				exit( EXIT_FAILURE );
			}

			std::string line;
			getline( in, line );

			TokType tok( line, boost::char_separator< char >( sep.c_str() ) );

			unsigned int results = 0;
			try
			{
				for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
					results++;
			}
			catch ( boost::bad_lexical_cast& e )
			{
				std::cout << "*** WARNING ***: bad lexical cast during training data parsing!" << std::endl;
			}
			in.close();
			return results;
		}

		/**
		 * Return number of training data points.
		 */
		unsigned int GetNumberOfTrainingPoints( const std::string& trainingData )
		{
			std::ifstream in( trainingData.c_str() );

			if ( in.fail() )
			{
				std::cerr << "*** ERROR ***: could not read labels from: " << trainingData << std::endl;
				exit( EXIT_FAILURE );
			}
			std::string line;
			unsigned int result = 0;

			while ( getline( in, line ) )
			{
				result++;
			}

			in.close();
			return result;
		}

		/**
		 * Add column with initial intercept values (1.0).
		 */
		MatrixType AddIntercept( const MatrixType& M )
		{
			MatrixType result( M.rows() + 1, M.cols() );

			for( unsigned int i = 0; i < result.rows(); i++ )
			{
				for( unsigned int j = 0; j < result.cols(); j++ )
				{
					if( i == 0 )
						result( i, j ) = 1.0;
					else
						result( i, j ) = M( i - 1, j );
				}
			}
			return result;
		}

		/**
		 * Read data from text file.
		 */
		void ReadData( const std::string& sep, const std::string& inputFileName, MatrixType& covariates, VectorType& binaryResponse, PixelType percentage )
		{
			unsigned int dims = GetNumberOfTrainingDims( inputFileName );
			unsigned int points = GetNumberOfTrainingPoints( inputFileName );

			covariates = MatrixType( dims - 1, points, 0.0 );
			binaryResponse = VectorType( points, 0.0 );

			// open training data ...
			std::ifstream in( inputFileName.c_str() );
			if ( in.fail() )
			{
				std::cerr << "*** ERROR ***: could not read data from: " << inputFileName << std::endl;
				exit( EXIT_FAILURE );
			}

			// insert all values into matrix, except the last one (responses) which is inserted in
			// the respose vector ...
			std::string line;
			unsigned int rowIndex = 0;
			while ( getline( in, line ) )
			{
				TokType tok( line, boost::char_separator< char >( sep.c_str() ) );
				try
				{
					unsigned int colIndex = 0;
					for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
					{
						if ( colIndex == dims - 1 ) // last column ...
							binaryResponse( rowIndex ) = boost::lexical_cast< PixelType >( *id );
						else
							covariates( colIndex, rowIndex ) = boost::lexical_cast< PixelType >( *id );

						colIndex++;
					}
					rowIndex++;
				} catch ( boost::bad_lexical_cast& e )
				{
					std::cout << "*** WARNING ***: could not parse data! (is it in column format?)" << std::endl;
					e.what();
					exit( EXIT_FAILURE );
				}
			}
			in.close();

			if ( ( percentage > 0 ) && ( percentage < 100 ) )
			{
				ResizeTrainingData( covariates, binaryResponse, percentage );
			}
		}

		/**
		 * Reduce number of training-points to given percentage.
		 */

		void ResizeTrainingData( MatrixType& trainingData, VectorType& trainingClasses,
				PixelType percentage )
		{
			unsigned int pointsToKeep = ( trainingClasses.size() * percentage ) / 100.0;

			// initialize vector ...
			VectorType subClasses( pointsToKeep, 0.0 );

			// initialize matrix ...
			MatrixType subData( trainingData.rows(), pointsToKeep, 0.0 );

			// create generator ...
			random_number_type generator( time( 0 ) );

			VectorType indices = GetRandomValues( trainingClasses, generator, pointsToKeep );

			// insert random samples ...
			for ( unsigned int i = 0; i < indices.size(); i++ )
			{
				subClasses( i ) = trainingClasses( indices( i ) );

				for ( unsigned int j = 0; j < subData.rows(); j++ )
					subData( j, i ) = trainingData( j, indices( i ) );
			}

			trainingData = subData;
			trainingClasses = subClasses;
		}

		/**
		 * Return uniform random values vector from inputs, with similar size.
		 */
		VectorType GetRandomValues( const VectorType& inputs, random_number_type& generator, unsigned int total )
		{
			// uniform distribution: 0 -> largest index ...
			int_distribution_type int_uni_dist( 0, inputs.size() - 1 );

			// generator ...
			int_generator_type int_distribution( generator, int_uni_dist );

			VectorType results( total, 0.0 );

			for ( unsigned int i = 0; i < total; i++ )
				results( int_distribution() );

			return results;
		}


		/**
		 * Return weight matrix.
		 */
		MatrixType Weights( const VectorType& means, const VectorType& coefficients, const MatrixType& covariate )
		{
			MatrixType weights( means.size(), means.size(), 0 );

			for ( unsigned int i = 0; i < means.size(); i++ )
				weights( i, i ) = means( i ) * ( 1. - ( means( i ) ) );

			return weights;
		}

		/**
		 * Return means vector.
		 */
		VectorType Means( const VectorType& coefficients, const MatrixType& covariate )
		{
			VectorType linearPredictors = GetColumnPackedCopy( covariate.transpose() * Vector2Matrix( coefficients ) );

			VectorType means( linearPredictors.size() );

			for ( unsigned int i = 0; i < means.size(); i++ )
				means( i) = std::exp( linearPredictors( i ) ) / ( static_cast< PixelType > ( 1.0 ) + std::exp( linearPredictors( i ) ) );

			return means;
		}

		/**
		 * Construct matrix from vector.
		 */
		MatrixType Vector2Matrix( const VectorType& V )
		{
			MatrixType A( V.size(), 1 );
			A.set_column( 0, V );

			return A;
		}

		/**
		 * Convert all matrix columns to vector.
		 */
		VectorType GetColumnPackedCopy( const MatrixType& M )
		{
			VectorType V( M.rows() * M.cols() );

			for ( unsigned int i = 0; i < M.rows(); i++ )
				for ( unsigned int j = 0; j < M.cols(); j++ )
					V( i + j * M.rows() ) = M( i, j );

			return V;
		}

		/**
		 * Set initial values for IRLS estimation.
		 */
		VectorType SetInit( const MatrixType& covariate )
		{
			VectorType V( covariate.rows() );

			PixelType a = 0;

			for ( unsigned int i = 0; i < V.size(); i++ )
			{
				a = 0;
				for ( unsigned int j = 0; j < covariate.get_row( i ).size(); j++ )
				{
					a += covariate( i, j );
				}

				V( i) = covariate.get_row( i ).size() / ( static_cast< PixelType > ( 100 ) * a );
			}
			return V;
		}

		/**
		 * Reweighted Least Squares ( IRLS ) estimation, by finding the
		 * maximum likelihood estimates of a generalized linear model( GLM ).
		 */
		VectorType IRLS( const VectorType& response, const MatrixType& covariate, unsigned int max_iter, bool verbose )
		{
			if ( response.size() != covariate.cols() )
			{
				std::cerr << "*** ERROR ***: The response vector and rows of the "
					"covariate matrix must have the same length." << std::endl;
				exit( EXIT_FAILURE );
			}

			VectorType coefficients = SetInit( covariate );

			MatrixType responseMatrix = Vector2Matrix( response );

			PixelType error = 1.0;

			MatrixType covariateMatrix( covariate );

			if( verbose )
				std::cout << "Starting iterative solving..." << std::endl;

			unsigned int iter = 0;
			while ( error > 0.000001 )
			{
				if( iter > max_iter )
				{
					std::cerr << "*** ERROR ***: maximum iteration reached. No convergence."<< std::endl;
					exit( EXIT_FAILURE );
				}

				if( verbose )
					std::cout << "Iteration: " << iter << std::endl;

				MatrixType coefficientMatrix = Vector2Matrix( coefficients );

				VectorType means = Means( coefficients, covariate );
				MatrixType weights = Weights( means, coefficients, covariate );

				MatrixType inversedWeights( weights.rows(), weights.cols() );
				MatrixType weight = MatrixType( weights );

				for ( unsigned int i = 0; i < weights.rows(); i++ )
				{
					if ( weights( i, i ) == 0. )
						inversedWeights( i, i ) = 0.;
					else
						inversedWeights( i, i ) = static_cast< PixelType > ( 1. ) / weights( i, i );
				}

				MatrixType linearPredictor = covariateMatrix.transpose() * coefficientMatrix;
				MatrixType z = linearPredictor + ( inversedWeights ) * ( responseMatrix - ( Vector2Matrix( means ) ) );
				MatrixType xwx = covariateMatrix * ( weight * covariateMatrix.transpose() );

				// Correct for rows with all values close to zero, by adding 0.1 to diagonal ...
				if ( std::abs( vnl_determinant< PixelType > ( xwx ) ) <= 1e-8 )
					for ( unsigned int i = 0; i < xwx.rows(); i++ )
						xwx( i, i ) = xwx( i, i ) + 0.1;

				// inverse...
				vnl_svd< PixelType > svd( xwx );

				MatrixType updatedCoefficient = svd.inverse() * ( ( covariateMatrix * weight ) * z );

				coefficients = GetColumnPackedCopy( updatedCoefficient );
				error = std::pow( ( updatedCoefficient - coefficientMatrix ).frobenius_norm(), 2.0 );

				iter++;
			}

			return coefficients;
		}

	};
} // end namespace generalized_linear_model



/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Generalized linear model" );

	parameters list;

	list.iterations = 100000;
	list.link = "logistic";
	list.separation = "\t";
	list.percentage = 100.0;
	list.verbose = false;
	list.printPValues = false;
	list.printZValues = false;

	p.AddArgument( list.inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription(
			"Input file: list of covariate(s) and response function (column format)" ) ->SetRequired( true );

	p.AddArgument( list.link, "link" ) ->AddAlias( "l" ) ->SetInput( "string" ) ->SetDescription( "Link function (default: \"logistic\")" );

	p.AddArgument( list.separation, "separation" ) ->AddAlias( "s" ) ->SetInput( "string" ) ->SetDescription( "Data separator string (default: \\t" );

	p.AddArgument( list.printZValues, "print-z-values" ) ->AddAlias( "pz" ) ->SetInput( "bool" ) ->SetDescription( "Print Z-values for coefficients (default: false" );

	p.AddArgument( list.printPValues, "print-p-values" ) ->AddAlias( "pp" ) ->SetInput( "bool" ) ->SetDescription( "Print P-values for coefficients (default: false" );

	p.AddArgument( list.percentage, "percentage" ) ->AddAlias( "p" ) ->SetInput( "float" ) ->SetDescription(
			"Percentage of covariates to select randomly (default: 100)" );

	p.AddArgument( list.verbose, "verbose" ) ->AddAlias( "v" ) ->SetInput( "bool" ) ->SetDescription( "Verbose (default: false" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	generalized_linear_model::GLM glm;

	glm.Run( list );

	return EXIT_SUCCESS;
}

