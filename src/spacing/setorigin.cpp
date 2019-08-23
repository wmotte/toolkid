/*
 * setspacing.cpp
 *
 *  Created on: Jun 10, 2009
 *      Author: wim
 */

#include "itkImageRegionConstIterator.h"
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

/**
 * In this class origin is set.
 */
class SetOrigin {

public:

	/**
	 * Write output after resetting origins.
	 */
	void setOrigin( const std::string& inputImageName, std::vector<double>& origins ) {

		// get dimensions...
		unsigned int dims = getDimensions( inputImageName );

		if ( dims < 2 || dims > 4 ) {
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		if ( origins.size() != dims )
		{
			std::cerr << "*** ERROR *** Image dimension (" << dims << "), does not "
				"match number of origins (" << origins.size() << ")." << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( inputImageName, origins ); \
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

	/**
	 * Process input image.
	 */
	template< class TPixel, unsigned int VDimension >
	void process( const std::string& inputImageName, const std::vector<double>& origins ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		// get input image...
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputImageName );
		typename ImageType::Pointer input = reader -> GetOutput();
		reader -> Update();

		// get origin...
		typename ImageType::PointType origin = input -> GetOrigin();

		for ( unsigned int i=0; i<origins.size(); i++ ) {
			origin[i] = origins[i];
		}

		input -> SetOrigin( origin );

		typename WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( inputImageName );
		writer -> SetInput( input );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error updating: " << inputImageName << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}
	}
};

/**
 * Set voxel spacing with given values. Update Input.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::vector<double> origins;
	std::string inputImageName;

	tkd::CmdParser
			parser( argv[0], "Set origin. (supports 2D, 3D and 4D images)" );

	parser.AddArgument( inputImageName, "input" )
		-> AddAlias( "i" )
		-> SetInput( "<string>" )
		-> SetDescription( "Input Image File" )
		-> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( origins, "origins" )
		-> AddAlias( "o" )
		-> SetInput( "<doubles>" )
		-> SetDescription( "Origins (e.g. -20  -18.2  5.3)" )
		-> SetRequired( true ) -> SetMinMax( 2, 4 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	SetOrigin setOrigin = SetOrigin();

	setOrigin.setOrigin( inputImageName, origins );

}

