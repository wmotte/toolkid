/*
 *
 *  Created on: Jul 09, 2011
 *      Author: wim
 */

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "tkdCmdParser.h"

/**
 * Return image size.
 */
class ImageSize {

public:

	/**
	 * Run.
	 */
	void run( const std::string& input ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( input ); \
		}

		switchMacro( float, 2 );
		switchMacro( float, 3 );
		switchMacro( float, 4 );

		exit( EXIT_SUCCESS );
	}

	/**
	 * Check dimensions of inputfile. In case of error, EXIT_FAILURE is returned.
	 */
	unsigned int getDimensions( const std::string& inputFileName ) {
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io ) {
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			return EXIT_FAILURE;
		} else {
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	template< class TPixel, unsigned int VDimension>
	void process( const std::string& inputFileName ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;

		// read input...
		typename ReaderType::Pointer inputReader = ReaderType::New();
		inputReader-> SetFileName( inputFileName );
		inputReader-> Update();
        std::cout << inputReader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;

	tkd::CmdParser parser( argv[0], "Get image dimensions [supports 2D,3D and 4D]" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ImageSize app = ImageSize();

	app.run( input );
}

