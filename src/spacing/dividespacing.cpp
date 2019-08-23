#include "itkImageRegionConstIterator.h"
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

/**
 * In this class voxel size of one specific direction is divided by 2. (used to BET ge3d).
 */
class DivideSpacing {

public:

	/**
	 * Write output after dividing one specific voxel spacing with division factor.
	 */
	void divide( const std::string& inputImageName, double inputDivideValue, unsigned int direction ) {

		// get dimensions...
		unsigned int dims = getDimensions( inputImageName );

		if ( dims < 2 || dims > 4 ) {
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		if ( direction < 0 || direction >= dims )
		{
			std::cerr << "*** ERROR *** Direction "  << direction << " is not supported for input image." << std::endl;
			exit( EXIT_FAILURE );
		}


		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( inputImageName, inputDivideValue, direction ); \
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
	void process( const std::string& inputImageName, double inputDivideValue, int direction ) {

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

		if ( inputDivideValue > 0 ) {
			spacing[direction] = spacing[direction] / inputDivideValue;
			origin[direction] = origin[direction] / inputDivideValue;
		} else {
			std::cout << "*** WARNING *** Division Value is equal or lower than zero: (" << inputDivideValue << ")!" << std::endl;
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
	int direction;
	double inputDivideValue;
	std::string inputImageName;

	tkd::CmdParser
			parser( argv[0], "Divide voxel spacing with constant in one direction and update origin. (supports 2D, 3D and 4D images)" );

	parser.AddArgument( inputDivideValue, "division-value" ) -> AddAlias( "m" ) -> SetInput( "<double>" ) -> SetDescription(
			"Division Value" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( inputImageName, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( direction, "direction" ) -> AddAlias( "d" ) -> SetInput( "<int>" ) -> SetDescription(
			"Direction to divide voxel spacing; where (x=0,y=1,z=2,t=3)" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	DivideSpacing divideSpacing = DivideSpacing();

	divideSpacing.divide( inputImageName, inputDivideValue, (unsigned) direction );

}

