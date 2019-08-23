#include "itkImageRegionConstIterator.h"
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

/**
 * In this class voxel sizes could be multiplied with a constant (e.g. 0.1 or 10.0).
 */
class ChangeSpacing {

public:

	/**
	 * Write output after multiplying voxel spacing with multiplication factor.
	 */
	void multiply( const std::string& inputImageName, double inputMultiplyValue ) {

		// get dimensions...
		unsigned int dims = getDimensions( inputImageName );

		if ( dims < 2 || dims > 4 ) {
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( inputImageName, inputMultiplyValue ); \
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
	void process( const std::string& inputImageName, double inputMultiplyValue ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		// get input image...
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputImageName );
		typename ImageType::Pointer input = reader -> GetOutput();
		reader -> Update();

		// get input spacing...
		typename ImageType::SpacingType spacing = input -> GetSpacing();

		// get origin...
		typename ImageType::PointType origin = input -> GetOrigin();

		// loop over spacing values...
		for ( unsigned int i = 0; i < VDimension; i++ )
		{
			double space = spacing[i] * inputMultiplyValue;
			spacing[i] = space;
		}

		// loop over origin...
		for ( unsigned int i = 0; i < VDimension; i++ )
		{
			double orig = origin[i] * inputMultiplyValue;
			origin[i] = orig;
		}

		input -> SetSpacing( spacing );
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
 * Multiply voxel spacing with given value. Update Input.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	double inputMultiplyValue;
	std::string inputImageName;

	tkd::CmdParser parser( argv[0], "Multiply voxel spacing with constant and update origin. (supports 2D, 3D and 4D images)" );

	parser.AddArgument( inputMultiplyValue, "multiply-value" ) -> AddAlias( "m" ) -> SetInput( "<double>" ) -> SetDescription(
			"Multiply Value" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( inputImageName, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ChangeSpacing changeSpacing = ChangeSpacing();

	changeSpacing.multiply( inputImageName, inputMultiplyValue );

}

