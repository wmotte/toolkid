#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkConnectedComponentImageFilter.h"

/**
 * Relabel binary Nissl stained images ...
 */
class Relabel
{

public:

	void run( const std::string& input, const std::string& output )
	{
		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 2 )
		{
			process2D( input, output );
		} else
		{
			std::cerr << "Number of dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Check dimensions of inputfile. In case of error, Application is terminated.
	 */
	unsigned int getDimensions( const std::string& inputFileName )
	{
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io )
		{
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else
		{
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Process 2D input image.
	 */
	void process2D( const std::string& input, const std::string& output )
	{
		typedef itk::Image< unsigned int, 2 > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		typedef itk::ConnectedComponentImageFilter< ImageType, ImageType > ConnectedComponentImageFilterType;

		// read...
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( input );
		reader->Update();

		ConnectedComponentImageFilterType::Pointer filter = ConnectedComponentImageFilterType::New();

		filter->SetInput( reader->GetOutput() );
		filter->Update();

		// write...
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( output );
		writer->SetInput( filter->GetOutput() );
		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::string input;
	std::string output;

	tkd::CmdParser parser( argv[0], "Relabel binary 2D image" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "2D Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "2D Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Relabel relabel = Relabel();

	relabel.run( input, output );
}
