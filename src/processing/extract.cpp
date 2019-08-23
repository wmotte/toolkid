/*
 * extract.cpp
 *
 *  Created on: Jun 26, 2009
 *      Author: wim
 */

/*
 * convert.cpp
 *
 *  Created on: Jun 26, 2009
 *      Author: wim
 */
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"

/**
 * Convert input to output format. (e.g. mhd to nii.gz).
 */
class Extract {

public:

	void run( const std::string& input, const std::string& output, unsigned int volume, unsigned int range ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 4 )
		{
			process4D( input, output, volume, range );
		} else
		{
			std::cerr << "Number of dimensions should be 4!" << std::endl;
			exit( EXIT_FAILURE );
		}
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
	 * Process 4D input image.
	 */
	void process4D( const std::string& input, const std::string& output, unsigned int volume, unsigned int range ) {

		typedef itk::Image< float, 3 > Image3DType;
		typedef itk::Image< float, 4 > Image4DType;
		typedef itk::ImageFileReader< Image4DType > ReaderType;
		typedef itk::ImageFileWriter< Image3DType > Writer3DType;
		typedef itk::ExtractImageFilter< Image4DType, Image3DType > ExtractImageFilter3DType;
		typedef itk::ImageFileWriter< Image4DType > Writer4DType;
		typedef itk::ExtractImageFilter< Image4DType, Image4DType > ExtractImageFilter4DType;

		// read...
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		reader -> Update();

		// region...
		Image4DType::RegionType region = reader -> GetOutput() -> GetLargestPossibleRegion();

		Image4DType::IndexType index = region.GetIndex();
		Image4DType::SizeType size = region.GetSize();

		// set index...
		if ( volume < size[3] ) {
			index[3] = volume;
		} else {
			std::cerr << "Index out of bounds! ( " << volume << " > " << ( size[3] - 1 ) << " )" << std::endl;
			exit( EXIT_FAILURE );
		}

		if ( range == 1 ) {
			// set size to 1...
			size[3] = 0;
			region.SetSize( size );
			region.SetIndex( index );
			ExtractImageFilter3DType::Pointer filter = ExtractImageFilter3DType::New();
			filter -> SetInput( reader -> GetOutput() );
			filter -> SetExtractionRegion( region );

			// write...
			Writer3DType::Pointer writer = Writer3DType::New();
			writer -> SetFileName( output );
			writer -> SetInput( filter -> GetOutput() );
			try {
				writer -> Update();
			} catch ( itk::ExceptionObject& e ) {
				std::cerr << "Error writing: " << output << std::endl;
				std::cerr << e.GetDescription() << std::endl;
			}

		} else {
			// set size to range...
			size[3] = range;
			region.SetSize( size );
			region.SetIndex( index );
			ExtractImageFilter4DType::Pointer filter = ExtractImageFilter4DType::New();
			filter -> SetInput( reader -> GetOutput() );
			filter -> SetExtractionRegion( region );

			// write...
			Writer4DType::Pointer writer = Writer4DType::New();
			writer -> SetFileName( output );
			writer -> SetInput( filter -> GetOutput() );
			try {
				writer -> Update();
			} catch ( itk::ExceptionObject& e ) {
				std::cerr << "Error writing: " << output << std::endl;
				std::cerr << e.GetDescription() << std::endl;
			}
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
	int volume = 0;
	int range  = 1;

	tkd::CmdParser parser( argv[0], "Extract volume from time-series or slice from volume" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( volume, "volume" ) -> AddAlias( "v" ) -> SetInput( "<uint>" ) -> SetDescription(
			"Volume or Slice to extract (default: 0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( range, "range" ) -> AddAlias( "r" ) -> SetInput( "<uint>" ) -> SetDescription(
			"Range to extract (default: 1)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Extract extract = Extract();

	extract.run( input, output, (unsigned) volume, (unsigned) range );
}
