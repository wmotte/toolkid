/*
 * crop.cpp
 *
 *  Created on: Jun 22, 2009
 *      Author: wim
 */

#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

/**
 * In this class the input image will be cropped automatically, between all voxels
 * with value < minValue.
 */
class Crop {

public:

	/**
	 * Run.
	 */
	void run( const std::string& input, const std::string& output, double minValue ) {

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
			process< pixel, dimension >( input, output, minValue ); \
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
	void process( const std::string& input, const std::string& output, float minValue ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;
		typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractImageFilterType;

		// get input image...
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		typename ImageType::ConstPointer image = reader -> GetOutput();
		reader -> Update();

		// find index and size to crop...
		IteratorType it( image, image -> GetLargestPossibleRegion() );

		typename ImageType::IndexType startCropIndex = image -> GetLargestPossibleRegion().GetIndex();
		typename ImageType::IndexType endCropIndex = image -> GetLargestPossibleRegion().GetIndex();
		typename ImageType::SizeType outputSize = image -> GetLargestPossibleRegion().GetSize();

		// set start crop index to maximum...
		for ( unsigned int i = 0; i < VDimension; i++ ) {
			startCropIndex[i] = outputSize[i] - 1;
		}

		// loop over voxels. If voxel value > minimum value, update start and end crop indices...
		for ( ; !it.IsAtEnd(); ++it ) {

			typename ImageType::IndexType index = it.GetIndex();
			TPixel intensity = image -> GetPixel( index );

			if ( intensity > minValue ) {
				for ( unsigned int i = 0; i < VDimension; i++ ) {
					if ( endCropIndex[i] < index[i] ) {
						endCropIndex[i] = index[i];
					} else if ( startCropIndex[i] > index[i] ) {
						startCropIndex[i] = index[i];
					}
				}
			}
		}

		// set output size...
		for ( unsigned int i = 0; i < VDimension; i++ ) {
			outputSize[i] = ( endCropIndex[i] - startCropIndex[i] ) + 1;
		}

		// set crop region...
		typename ImageType::RegionType outputRegion;
		outputRegion.SetSize( outputSize );
		outputRegion.SetIndex( startCropIndex );

		// filter...
		typename ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
		extractFilter -> SetInput( image );
		extractFilter -> SetExtractionRegion( outputRegion );

		// write output...
		typename WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( output );
		writer -> SetInput( extractFilter -> GetOutput() );

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
	double minValue = 0.0;

	tkd::CmdParser parser( argv[0], "Crop input image. (supports 2D, 3D and 4D images)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( minValue, "min-value" ) -> AddAlias( "m" ) -> SetInput( "<double>" ) -> SetDescription(
			"Minimum value to crop image with (default: 0.0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Crop crop = Crop();

	crop.run( input, output, minValue );
}

