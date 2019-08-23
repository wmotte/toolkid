/*
 * multicrop.cpp
 *
 * Author: wim
 */

#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

/**
 * In this class the input images will be cropped automatically, between all voxels
 * with value < minValue.
 */
class MultiCrop {

public:

	/**
	 * Run.
	 */
	void run( const std::string& initial, const std::vector< std::string >& inputs, const std::vector< std::string >& outputs,
			double minValue ) {

		// check correct sizes!
		if ( inputs.size() != outputs.size() ) {
			std::cerr << "Inputs and Outputs should have same size!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// get dimensions...
		unsigned int dims = getDimensions( initial );

		if ( dims < 2 || dims > 4 ) {
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( initial, inputs, outputs, minValue ); \
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
	void process( const std::string& initial, const std::vector< std::string >& inputs,
					const std::vector< std::string >& outputs, float minValue )
	{
		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;
		typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractImageFilterType;

		// get input image...
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( initial );
		typename ImageType::ConstPointer image = reader -> GetOutput();
		reader -> Update();

		// find index and size to crop...
		IteratorType it( image, image -> GetLargestPossibleRegion() );

		typename ImageType::IndexType startCropIndex = image -> GetLargestPossibleRegion().GetIndex();
		typename ImageType::IndexType endCropIndex   = image -> GetLargestPossibleRegion().GetIndex();
		typename ImageType::SizeType outputSize      = image -> GetLargestPossibleRegion().GetSize();

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

		// crop all inputs...
		for ( unsigned int i = 0; i < inputs.size(); i++ ) {

			std::string input = inputs[i];
			std::string output = outputs[i];

			typename ReaderType::Pointer reader = ReaderType::New();
			reader -> SetFileName( input );
			typename ImageType::ConstPointer image = reader -> GetOutput();
			reader -> Update();

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
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string initial;
	std::vector< std::string > inputs;
	std::vector< std::string > outputs;

	double minValue = 0.0;

	tkd::CmdParser parser( argv[0], "Multi Crop images. (2D, 3D and 4D)" );

	parser.AddArgument( initial, "initial" ) -> AddAlias( "t" ) -> SetInput( "<string>" ) -> SetDescription(
			"Initial Image used to determine crop region" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( inputs, "inputs" ) -> AddAlias( "i" ) -> SetInput( "vector<string>" ) -> SetDescription( "Images" ) -> SetRequired(
			true );

	parser.AddArgument( outputs, "outputs" ) -> AddAlias( "o" ) -> SetInput( "vector<string>" ) -> SetDescription( "Outputs" ) -> SetRequired(
			true );

	parser.AddArgument( minValue, "min-value" ) -> AddAlias( "m" ) -> SetInput( "<double>" ) -> SetDescription(
			"Minimum value to crop initial image with (default: 0.0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	MultiCrop multicrop = MultiCrop();

	multicrop.run( initial, inputs, outputs, minValue );
}

