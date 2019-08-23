/*
 * indices.cpp
 *
 *  Created on: Sep 23, 2009
 *      Author: wim
 */
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

/**
 * Print row of all indices to text file of non-zero voxels.
 */
class Indices {

public:

	/**
	 * Run indices.
	 */
	void run( const std::string& input, const std::string& output ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 3 ) {
			process3D( input, output );
		} else {
			std::cout << "Dims != 3 not supported!" << std::endl;
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
	 * Process 3D input image.
	 */
	void process3D( const std::string& input, const std::string& output ) {

		typedef itk::Image< unsigned char, 3 > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

		// get input image...
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		ImageType::Pointer image = reader -> GetOutput();
		reader -> Update();

		ImageType::RegionType region = image -> GetLargestPossibleRegion();

		// find index and size to crop...
		IteratorType it( image, region );

		// write output...
		std::ofstream filestream;
		filestream.open( output.c_str() );
		if ( filestream.fail() ) {
			std::cerr << "Not able to write to: " << output << std::endl;
			exit( EXIT_FAILURE );
		}

		// set all voxels in clearRegion
		for ( ; !it.IsAtEnd(); ++it ) {

			ImageType::IndexType index = it.GetIndex();
			ImageType::PixelType pixel = image -> GetPixel( index );

			if ( pixel != 0 ) {
				filestream << index[0] << "," << index[1] << "," << index[2] << std::endl;
			}
		}

		filestream.close();
	}
}
;

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;
	std::string output;

	tkd::CmdParser parser( argv[0], "Print all indices from non-zero voxels" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Text File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Indices indices = Indices();

	indices.run( input, output );
}
