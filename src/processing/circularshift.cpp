#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "tkdCmdParser.h"

/**
 * In this class the input images is circular shifted.
 */
class CircularShift {

public:

	/**
	 * Run circular shift.
	 */
	void run( const std::string& input, const std::string& output, unsigned int direction, int shift ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		if( dims < 2 || dims > 4 ) {
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		if ( direction < 0 || direction >= dims )
		{
			std::cerr << "*** ERROR *** Direction is not supported for input image: " << direction << std::endl;
			exit( EXIT_FAILURE );
		}

		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( input, output, direction, shift ); \
		}

		switchMacro( float, 2 );
		switchMacro( float, 3 );
		switchMacro( float, 4 );
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
	 * Circular shift.
	 */
	template< class TPixel, unsigned int VDimension >
	void process( const std::string& inputFileName, const std::string& outputFileName, int direction, int shift )
	{
		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputFileName );
		reader -> Update();

		typename ImageType::Pointer image = reader->GetOutput();
		typename ImageType::Pointer output = ImageType::New();
		output -> CopyInformation( image );
		output -> SetRegions( image -> GetLargestPossibleRegion() );
		output -> Allocate();

		typename ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
		int dimension = size[direction];

		typename itk::ImageRegionIteratorWithIndex< ImageType > it( output, output->GetLargestPossibleRegion() );
		for ( ; !it.IsAtEnd(); ++it ) {
			typename ImageType::IndexType index = it.GetIndex();
			index[direction] = ( index[direction] + shift ) % dimension;

			it.Set( image -> GetPixel( index ) );
		}

		typename WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( outputFileName );
		writer -> SetInput( output );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << outputFileName << std::endl;
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
	int direction;
	int shift;

	tkd::CmdParser parser( argv[0], "Circular shift image in one direction given the shift in voxels. (supports 2D, 3D and 4D images)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( direction, "direction" ) -> AddAlias( "d" ) -> SetInput( "<int>" ) -> SetDescription(
			"Direction to shift ( x=0; y=1; z=2; t=3 )" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( shift, "shift" ) -> AddAlias( "s" ) -> SetInput( "<int>" ) -> SetDescription( "Number of voxels to shift" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	CircularShift circularShift = CircularShift();

	circularShift.run( input, output, (unsigned)direction, shift );
}

