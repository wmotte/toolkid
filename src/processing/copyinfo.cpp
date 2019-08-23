/*
 * copyinfo.cpp
 *
 *  Created on: Jul 21, 2009
 *      Author: wim
 */
#include "itkImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "tkdCmdParser.h"

#include <iostream>
#include <fstream>
#include <string>

/**
 * Copy header info from one file to another.
 */
class CopyInfo {

public:

	void run( const std::string& input, const std::string& reference ) {

		// get dimensions...
		unsigned int dimsInput = getDimensions( input );
		unsigned int dimsReference = getDimensions( reference );

		// inputs should have similar size...
		if ( dimsInput != dimsReference )
		{
			std::cerr << "Number of dimensions differs between input and reference!" << std::endl;
			exit( EXIT_FAILURE );
		}
		else if( dimsInput < 2 || dimsInput >4 )
		{
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dimsInput == dimension ) \
		{ \
			process< pixel, dimension >( input, reference ); \
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
	void process( const std::string& inputFileName, const std::string& referenceFileName ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageFileReader< ImageType > ReaderType;

		typename ReaderType::Pointer inputReader = ReaderType::New();
		typename ReaderType::Pointer referenceReader = ReaderType::New();
		typename WriterType::Pointer writer = WriterType::New();

		inputReader -> SetFileName( inputFileName );
		referenceReader -> SetFileName( referenceFileName );

		typename ImageType::Pointer input = inputReader -> GetOutput();
		typename ImageType::Pointer reference = referenceReader -> GetOutput();

		inputReader -> Update();
		referenceReader -> Update();

		input -> CopyInformation( reference );
		input -> SetOrigin( reference -> GetOrigin() );

		writer -> SetFileName( inputFileName );
		writer -> SetInput( input );

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
	std::string reference;
	std::string input;

	tkd::CmdParser parser( argv[0], "Copy header info from one file to another (2D, 3D, 4D supported)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( reference, "reference" ) -> AddAlias( "r" ) -> SetInput( "<string>" ) -> SetDescription( "Reference File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	CopyInfo copyInfo = CopyInfo();
	copyInfo.run( input, reference );
}

