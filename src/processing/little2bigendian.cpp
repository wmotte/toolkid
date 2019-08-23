#include "itkImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkBigMetaImageIO.h"
#include "tkdCmdParser.h"

/**
 * Swap from little to big endianess, both 3D and 4D.
 */
class Little2BigEndian {

public:

	void run( const std::string& input, const std::string& output ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );


		if ( dims < 2 || dims > 4 ) {
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( input, output ); \
		}

		switchMacro( float, 2 );
		switchMacro( float, 3 );
		switchMacro( float, 4 );

		exit( EXIT_SUCCESS );
	}

	/**
	 * Check dimensions of inputfile. In case of error, Application is terminated.
	 */
	unsigned int getDimensions( const std::string& inputFileName ) {
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io ) {
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else {
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Process input image.
	 */
	template< class TPixel, unsigned int VDimension >
	void process( const std::string& input, const std::string& output ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::BigMetaImageIO BigMetaImageType;

		typename BigMetaImageType::Pointer meta = BigMetaImageType::New();
		typename ReaderType::Pointer reader = ReaderType::New();
		typename WriterType::Pointer writer = WriterType::New();

		reader -> SetFileName( input );

		writer -> SetImageIO( meta );
		writer -> SetFileName( output );
		writer -> SetInput( reader -> GetOutput() );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject & exp ) {
			std::cout << "Error writing image: " << input << std::endl;
			exit( EXIT_FAILURE );
		}
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;
	std::string output;

	tkd::CmdParser parser( argv[0], "Swap little-endian images to big-endian Meta Images (2D, 3D and 4D supported)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Little2BigEndian little2BigEndian = Little2BigEndian();

	little2BigEndian.run( input, output );
}

