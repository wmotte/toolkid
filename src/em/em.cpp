#include "tkdCmdParser.h"

#include "graphCommon.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSTAPLEImageFilter.h"


/**
 * wim@invivonmr.uu.nl, 14-08-2010.
 */
namespace em
{
	typedef double PixelType;
	typedef std::vector< PixelType > VectorType;

	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
	typedef itk::ImageFileWriter< ImageType > ImageFileWriterType;
	typedef itk::STAPLEImageFilter< ImageType, ImageType > EMFilterType;

	// *******************************************************************************

	/**
	 * Option list.
	 */
	struct parameters
	{
		std::vector< std::string > inputFileNames;
		std::string outputFileName;
		bool verbose;
	};

	// *******************************************************************************

	/**
	 *
	 */
	class EM
	{

	public:

		/**
		 * Run GLMM select voxel application.
		 */
		void Run( parameters& args )
		{
			Checks( args );

			EMFilterType::Pointer em = EMFilterType::New();

			for( unsigned int i = 0; i < args.inputFileNames.size(); i++ )
				em->SetInput( i, ReadImage( args.inputFileNames.at( i ) ) );


			em->Update();

			WriteImage( em->GetOutput(), args.outputFileName );

			if( args.verbose )
			{
				for( unsigned int i = 0; i < args.inputFileNames.size(); i++ )
				{
					std::cout << "Sensitivity: " << args.inputFileNames.at( i ) << " -> " << em->GetSensitivity( i ) << std::endl;
					std::cout << "Specificity: " << args.inputFileNames.at( i ) << " -> " << em->GetSpecificity( i ) << std::endl;
					std::cout << std::endl;
				}
				std::cout << "Confidence weight: " << em->GetConfidenceWeight() << std::endl;
				std::cout << "Elapsed iterations: " << em->GetElapsedIterations() << std::endl;

			}
		}

	private:

		/**
		 * Write.
		 */
		void WriteImage( const ImageType::Pointer image, const std::string outputFileName )
		{
			ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
			writer->SetFileName( outputFileName );
			writer->SetInput( image );

			try
			{
				writer->Update();
			} catch( ... )
			{
				std::cout << "*** ERROR ***: unable to write to " << outputFileName << std::endl;
				exit( EXIT_FAILURE );
			}
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
		 * Checks.
		 */
		void Checks( parameters& args )
		{
			// dims
			for( unsigned int i = 0; i < args.inputFileNames.size(); i++ )
			{
				if( graph::Graph< PixelType >::GetImageDimensions( args.inputFileNames.at( i ) ) != 3 )
				{
					std::cerr << "*** ERROR ***: not all images are 3D!" << std::endl;
					exit( EXIT_FAILURE );
				}
			}
		}

	};

} // end namespace em

/**
 * Main; expectation maximization.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Expectation maximization binary segmentations." );

	em::parameters args;
	args.verbose = false;

	p.AddArgument( args.inputFileNames, "inputs" )
			->AddAlias( "i" )
			->SetDescription( "3D Segmentation input files" )
			->SetMinMax( 1,10000 )
			->SetRequired( true );

	p.AddArgument( args.outputFileName, "output" )
			->AddAlias( "o" )
			->SetDescription( "3D Expectation-Maximization output file" )
			->SetMinMax( 1,1 )
			->SetRequired( true );

	p.AddArgument( args.verbose, "verbose" )
			->AddAlias( "v" )
			->SetDescription( "Verbose" );



	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	em::EM em;

	em.Run( args );

	return EXIT_SUCCESS;
}


