/*
 * rescale.cpp
 *
 *  Created on: Seo 26, 2009
 *      Author: wim
 */

#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

/**
 * Rescale input intensities between min and max values.
 */
class Rescale {

public:

	void run( const std::string& input, const std::string& output, float min, float max ) {

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
			process< pixel, dimension >( input, output, min, max ); \
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
	void process( const std::string& input, const std::string& output, float min, float max ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleIntensityImageFilterType;

		// read...
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		reader -> Update();

		typename RescaleIntensityImageFilterType::Pointer filter = RescaleIntensityImageFilterType::New();
		filter -> SetInput( reader -> GetOutput() );
		filter -> SetOutputMinimum( min );
		filter -> SetOutputMaximum( max );

		// write...
		typename WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( output );
		writer -> SetInput( filter -> GetOutput() );
		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
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
	float min = 0;
	float max = 100;

	tkd::CmdParser parser( argv[0], "Rescale image intensities between min and max (default 0 - 100)." );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( min, "min" ) -> AddAlias( "l" ) -> SetInput( "<float>" ) -> SetDescription( "Minimum of range" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( max, "max" ) -> AddAlias( "u" ) -> SetInput( "<float>" ) -> SetDescription( "Maximum of range" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Rescale rescale = Rescale();

	rescale.run( input, output, min, max );
}
