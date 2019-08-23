/*
 * checkerboard.cpp
 *
 *  Created on: Jun 11, 2009
 *      Author: wim
 *  Added templates: 11-11-2009
 */

#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkNormalizeImageFilter.h"

/**
 * In this class the checkerboard with normalize filter is used to craete checkerboards.
 */
class CheckerBoard {

public:

	/**
	 * Write checkerboard.
	 */
	void run( const std::vector< std::string >& inputs, const std::string& output, const std::vector< int > pattern )
	{
		// get dimensions...
		unsigned int dims1 = getDimensions( inputs[0] );
		unsigned int dims2 = getDimensions( inputs[1] );

		if( dims1 != dims2 ) {
			std::cerr << "Dimensions of input files are different! " << dims1 << "!=" << dims2 << std::endl;
			exit( EXIT_FAILURE );
		}
		else if( dims1 != pattern.size() )
		{
			std::cerr << "Number of given patterns is not equal to dimensions of input files! "
			<< dims1 << "!=" << pattern.size() << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims1 == dimension ) \
		{ \
			process< pixel, dimension >( inputs, output, pattern ); \
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
	void process( const std::vector< std::string >& inputs, const std::string& output, const std::vector< int > pattern ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		typedef itk::NormalizeImageFilter< ImageType, ImageType > NormalizeImageFilterType;
		typedef itk::CheckerBoardImageFilter< ImageType > CheckerBoardImageFilterType;

		typedef itk::FixedArray< unsigned int, VDimension > PatternType;

		// get input images...
		typename ReaderType::Pointer reader1 = ReaderType::New();
		typename ReaderType::Pointer reader2 = ReaderType::New();

		reader1 -> SetFileName( inputs[0] );
		reader2 -> SetFileName( inputs[1] );

		// normalize inputs...
		typename NormalizeImageFilterType::Pointer normalize1 = NormalizeImageFilterType::New();
		typename NormalizeImageFilterType::Pointer normalize2 = NormalizeImageFilterType::New();

		normalize1 -> SetInput( reader1 -> GetOutput() );
		normalize2 -> SetInput( reader2 -> GetOutput() );

		// create checker...
		typename CheckerBoardImageFilterType::Pointer checker = CheckerBoardImageFilterType::New();

		checker -> SetInput1( normalize1 -> GetOutput() );
		checker -> SetInput2( normalize2 -> GetOutput() );

		PatternType itkPattern;
		for ( unsigned int i = 0; i < pattern.size(); i++ ) {
			itkPattern[i] = pattern[i];
		}

		checker -> SetCheckerPattern( itkPattern );
		checker -> Update();

		// write output...
		typename WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( output );
		writer -> SetInput( checker -> GetOutput() );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}
	}

	/**
	 * Process 3D input image.
	 */
	void process3D( const std::vector< std::string>& inputs, const std::string& output, const std::vector< int> pattern ) {

		const unsigned int Dimension = 3;
		typedef float PixelType;
		typedef itk::Image< PixelType, Dimension> ImageType;
		typedef itk::ImageFileReader< ImageType> ReaderType;
		typedef itk::ImageFileWriter< ImageType> WriterType;

		typedef itk::NormalizeImageFilter< ImageType, ImageType> NormalizeImageFilterType;
		typedef itk::CheckerBoardImageFilter< ImageType> CheckerBoardImageFilterType;

		// get input images...
		ReaderType::Pointer reader1 = ReaderType::New();
		ReaderType::Pointer reader2 = ReaderType::New();

		reader1 -> SetFileName( inputs[0] );
		reader2 -> SetFileName( inputs[1] );

		// normalize inputs...
		NormalizeImageFilterType::Pointer normalize1 = NormalizeImageFilterType::New();
		NormalizeImageFilterType::Pointer normalize2 = NormalizeImageFilterType::New();

		normalize1 -> SetInput( reader1 -> GetOutput() );
		normalize2 -> SetInput( reader2 -> GetOutput() );

		// create checker...
		CheckerBoardImageFilterType::Pointer checker = CheckerBoardImageFilterType::New();

		checker -> SetInput1( normalize1 -> GetOutput() );
		checker -> SetInput2( normalize2 -> GetOutput() );

		itk::FixedArray< unsigned int, Dimension> itkPattern;
		for ( unsigned int i = 0; i < pattern.size(); i++ ) {
			itkPattern[i] = pattern[i];
		}

		checker -> SetCheckerPattern( itkPattern );
		checker -> Update();

		// write output...
		WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( output );
		writer -> SetInput( checker -> GetOutput() );

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
	std::string output;
	std::vector< std::string > inputs;
	std::vector< int > pattern;

	tkd::CmdParser parser( argv[0], "Create checkerboard image. (supports 2D, 3D and 4D images)" );

	parser.AddArgument( inputs, "inputs" ) -> AddAlias( "i" ) -> SetInput( "<strings>" ) -> SetDescription( "Two Input Image Files" ) -> SetRequired(
			true ) -> SetMinMax( 2, 2 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( pattern, "pattern" ) -> AddAlias( "p" ) -> SetInput( "<ints>" ) -> SetDescription(
			"Number of checks to make per image dimension" ) -> SetRequired( true ) -> SetMinMax( 2, 4 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	CheckerBoard checkerBoard = CheckerBoard();

	checkerBoard.run( inputs, output, pattern );
}

