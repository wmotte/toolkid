#include "tkdCmdParser.h"


#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "graphCommon.h"

#include "itkImageIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>



/**
 * wim@invivonmr.uu.nl, 31-07-2010.
 */
namespace glmm
{
	typedef double PixelType;
	typedef std::vector< PixelType > VectorType;
	typedef std::vector< VectorType > DataContainerType;

	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;

	typedef boost::mt19937 random_number_type;
	typedef boost::uniform_int< int > int_distribution_type;
	typedef boost::variate_generator< random_number_type&, int_distribution_type> int_generator_type;



	// *******************************************************************************

	/**
	 * Option list.
	 */
	struct glmmParameters
	{
		std::vector< std::string > inputFileNames;
		std::vector< std::string > labelFileNames;
		std::string outputFileName;

		std::vector< std::string > inputNames;
		std::vector< std::string > labelNames;
		std::vector< std::string > covariateNames;
		std::vector< PixelType > covariateValues;

		int seed;

		bool allVoxels;
	};

	// *******************************************************************************

	/**
	 *
	 */
	class GLMMSelectVoxels
	{

	public:

		/**
		 * Run GLMM select voxel application.
		 */
		void Run( glmmParameters& args )
		{
			// Check inputs.
			Checks( args );

			// Read inputs.
			ReadImages( args.inputFileNames );

			// Read labels.
			ReadLabels( args.labelFileNames );

			// Create datacontainer.
			ReadData();

			// Get minimal vector size and index.
			std::vector< unsigned int > sizes;
			for( unsigned int label = 0; label < m_Data.size(); label++ )
				sizes.push_back( ( m_Data[ label ][ 0 ] ).size() );

			std::vector< unsigned int >::iterator min_it = std::min_element( sizes.begin(), sizes.end() );

			unsigned int minIndex = std::distance( sizes.begin(), min_it );
			unsigned int minSize = sizes[ minIndex ];

			// Resize all matrices so all have equal size.
			if( ! args.allVoxels )
			{
				ResampleDataContainer( minSize, minIndex, args.seed );
			}

			// Add optional covariates
			AddCovariates( args.covariateValues );

			// Write datacontainer to csv.
			WriteCSV( args );

			exit( EXIT_SUCCESS );
		}

	private:

		std::vector< ImageType::Pointer > m_Images;
		std::vector< ImageType::Pointer > m_Labels;
		std::vector< DataContainerType > m_Data;

		/**
		 * Add covariates to data containers.
		 */
		void AddCovariates( const std::vector< PixelType >& values )
		{
			for( unsigned int m = 0; m < m_Data.size(); m++ )
			{
				for( unsigned c = 0; c < values.size(); c++ )
				{
					VectorType covariate( ( m_Data.at( m ) ).at( 0 ).size(), values.at( c ) );
					m_Data.at( m ).push_back( covariate );
				}
			}
		}

		/**
		 * Write matrices to csv.
		 */
		void WriteCSV( const glmmParameters& args )
		{
			// open stream
			std::ofstream out( args.outputFileName.c_str() );
			if ( out.fail() )
			{
				std::cerr << "*** ERROR ***: Could not write to '" << args.outputFileName << "'!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// write header
			std::vector< std::string > names;
			names.push_back( "Label" );

			std::back_insert_iterator< std::vector< std::string > > p( names );

			copy( args.inputNames.begin( ), args.inputNames.end( ), p );
			copy( args.covariateNames.begin( ), args.covariateNames.end( ), p );

			for( unsigned int i = 0; i < names.size(); i++ )
			{
				out << names.at( i );
				if( i < names.size() - 1 )
					out << ",";
			}
			out << std::endl; // end header

			// write data

			for( unsigned l = 0; l < m_Data.size(); l++ )
			{
				DataContainerType data = m_Data.at( l );

				unsigned int cols = data.size();
				unsigned int rows = data.at( 0 ).size();

				for( unsigned r = 0; r < rows; r++ )
				{
					out << args.labelNames.at( l ) << ",";
					for( unsigned c = 0; c < cols; c++ )
					{
						out << data[c][r];

						if( c < cols - 1 )
							out << ",";
					}
					out << std::endl;
				}
			}

			out.close();
		}

		/**
		 * Resample matrices, so they match the smallest one in size.
		 */
		void ResampleDataContainer( unsigned int minSize, unsigned int minIndex, int seed )
		{
			for( unsigned int label = 0; label < m_Data.size(); label++ )
			{
				if( label == minIndex )
					continue; // skip subsampling smallest label
				else
				{
					DataContainerType data = m_Data.at( label );

					for( unsigned int img = 0; img < data.size(); img++ )
						data.at( img ) = GetRandomValues( data.at( img ), minSize, seed );

					m_Data.at( label ) = data;
				}
			}
		}

		/**
		 * Read data containers (for each image one row) for each label image.
		 */
		void ReadData()
		{
			for( unsigned int label = 0; label < m_Labels.size(); label++ )
			{
				DataContainerType data;

				for( unsigned int image = 0; image < m_Images.size(); image++ )
				{
					VectorType vector;

					IteratorType it( m_Images.at( image ), m_Images.at( image )->GetLargestPossibleRegion() );

					for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
					{
						PixelType takePixel = m_Labels.at( label )->GetPixel( it.GetIndex() );

						if ( takePixel != 0 )
							vector.push_back( it.Get() );
					}
					data.push_back( vector );
				}
				m_Data.push_back( data );
			}
		}

		/**
		 * Allocate and read global images.
		 */
		void ReadImages( const std::vector< std::string >& images )
		{
			for( unsigned int i = 0; i < images.size(); i++ )
				m_Images.push_back( ReadImage( images.at( i ) ) );
		}

		/**
		 * Allocate and read global labels.
		 */
		void ReadLabels( const std::vector< std::string >& images )
		{
			for( unsigned int i = 0; i < images.size(); i++ )
				m_Labels.push_back( ReadImage( images.at( i ) ) );
		}

		/**
		 * Read image.
		 */
		ImageType::Pointer ReadImage( const std::string& image )
		{
			ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
			reader->SetFileName( image );
			reader->Update();
			return reader->GetOutput();
		}

		/**
		 * Return uniform random values vector from inputs, with similar size.
		 */
		VectorType GetRandomValues( const VectorType& input, unsigned int total, int seed )
		{
			// create generator ...
			random_number_type generator( seed );

			VectorType result;

			// uniform distribution: 0 -> largest index ...
			int_distribution_type int_uni_dist( 0, input.size() - 1 );

			// generator ...
			int_generator_type int_distribution( generator, int_uni_dist );

			for ( unsigned int i = 0; i < total; i++ )
			{
				//result.push_back( int_distribution() );
				result.push_back( input.at( int_distribution() ) );
			}
			return result;
		}

		/**
		 * Checks.
		 */
		void Checks( glmmParameters& args )
		{
			// if names not specified, use fileNames instead.
			if( args.inputNames.empty() )
				args.inputNames = args.inputFileNames;

			if( args.labelNames.empty() )
				args.labelNames = args.labelFileNames;

			// input names size
			if( args.inputNames.size() != args.inputFileNames.size() )
			{
				std::cerr << "*** ERROR ***: input names and input files different is size!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// label names size
			if( args.labelNames.size() != args.labelFileNames.size() )
			{
				std::cerr << "*** ERROR ***: label names and label files different is size!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// covariate size
			if( args.covariateNames.size() != args.covariateValues.size() )
			{
				std::cerr << "*** ERROR ***: covariate names and values different is size!" << std::endl;
				exit( EXIT_FAILURE );
			}
			// image dims
			for( unsigned int i = 0; i < args.inputFileNames.size(); i++ )
			{
				if( graph::Graph< PixelType >::GetImageDimensions( args.inputFileNames.at( i ) ) != 3 )
				{
					std::cerr << "*** ERROR ***: not all input images are 3D!" << std::endl;
					exit( EXIT_FAILURE );
				}
			}
			// labels dims
			for( unsigned int i = 0; i < args.labelFileNames.size(); i++ )
			{
				if( graph::Graph< PixelType >::GetImageDimensions( args.labelFileNames.at( i ) ) != 3 )
				{
					std::cerr << "*** ERROR ***: not all label images are 3D!" << std::endl;
					exit( EXIT_FAILURE );
				}
			}
		}

	};

} // end namespace glmm

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Create csv output for GLMM R." );

	glmm::glmmParameters args;

	args.seed = 100;
	args.allVoxels = false;

	p.AddArgument( args.inputFileNames, "inputs" )
			->AddAlias( "i" )
			->SetDescription( "3D Input files" )
			->SetMinMax( 1,10000 )
			->SetRequired( true );

	p.AddArgument( args.labelFileNames, "labels" )
			->AddAlias( "l" )
			->SetDescription( "3D Label files" )
			->SetMinMax( 1, 10000 )
			->SetRequired( true );

	p.AddArgument( args.outputFileName, "output" )
			->AddAlias( "o" )
			->SetDescription( "output CSV file" )
			->SetMinMax( 1,1 )
			->SetRequired( true );

	p.AddArgument( args.inputNames, "input-names" )
			->AddAlias( "in" )
			->SetDescription( "3D Input file names (default: inputFileNames)" )
			->SetMinMax( 1,10000 );

	p.AddArgument( args.labelNames, "label-names" )
			->AddAlias( "ln" )
			->SetDescription( "3D Label file names (default: labelFileNames)" )
			->SetMinMax( 1,10000 );

	p.AddArgument( args.covariateNames, "covariate-names" )
			->AddAlias( "cn" )
			->SetDescription( "Names optional covariates" )
			->SetMinMax( 1,10000 );

	p.AddArgument( args.covariateValues, "covariate-values" )
			->AddAlias( "cv" )
			->SetDescription( "Values optional covariates" )
			->SetMinMax( 1,10000 );

	p.AddArgument( args.seed, "seed" )
			->AddAlias( "s" )
			->SetDescription( "Seed (default: 100" )
			->SetMinMax( 1,1 );

	p.AddArgument( args.allVoxels, "all-voxels" )
			->AddAlias( "a" )
			->SetDescription( "Include all voxels (default: equal number of voxels per label" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	glmm::GLMMSelectVoxels glmmSelectVoxels;

	glmmSelectVoxels.Run( args );

	return EXIT_SUCCESS;
}


