/*
 * mean.cpp
 *
 *  Created on: Jun 26, 2009
 *      Author: wim
 */

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "tkdCmdParser.h"
//#include <numeric>

/**
 * In this class mean is calculated between bounds.
 */
class Mean {

public:

	/**
	 * Run circular shift.
	 */
	void run( const std::string& input, const std::string& mask, float lower, float upper ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( input, mask, lower, upper ); \
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

	template< class TPixel, unsigned int VDimension>
	void process( const std::string& inputFileName, const std::string& maskFileName, float lower, float upper ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;

		// read input...
		typename ReaderType::Pointer inputReader = ReaderType::New();
		inputReader -> SetFileName( inputFileName );
		inputReader -> Update();

		// read mask...
		typename ReaderType::Pointer maskReader = ReaderType::New();
		maskReader -> SetFileName( maskFileName );
		maskReader -> Update();

		typename ImageType::ConstPointer input = inputReader -> GetOutput();
		typename ImageType::ConstPointer mask  = maskReader  -> GetOutput();

		typename ImageType::RegionType inputRegion = input -> GetLargestPossibleRegion();
		typename ImageType::RegionType maskRegion  = mask  -> GetLargestPossibleRegion();

		typename ImageType::SizeType inputSize = inputRegion.GetSize();
		typename ImageType::SizeType maskSize  = maskRegion.GetSize();

		// size check...
		for ( unsigned int i=0; i<VDimension; i++ ) {
			if ( inputSize[i] != maskSize[i] ) {
				std::cerr << "Image and Mask Size don't match!" << inputSize << " != " << maskSize << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		std::vector< TPixel > includedValues;
		IteratorType it( mask, maskRegion );

		// iterate over mask; if mask != 0, inspect image value.
		// if lower <= value <= upper, include value in includedValues vector.
		for ( ; !it.IsAtEnd(); ++it ) {
			typename ImageType::IndexType index = it.GetIndex();
			TPixel maskPixel = mask -> GetPixel( index );

			if ( maskPixel != 0 ) {
				TPixel inputPixel = input -> GetPixel( index );

				if ( ( inputPixel >= lower ) && ( inputPixel <= upper ) ) {
					includedValues.push_back( inputPixel );
				}
			}
		}

		// output...
		if ( ! includedValues.empty() ) {

			TPixel sum = 0;
			for ( unsigned int i=0; i< includedValues.size(); i++ ) {
				sum += includedValues[i];
			}

			float average = ( sum / includedValues.size() );
			std::cout << average << std::endl;

		} else {
      std::cout << 0 << std::endl;
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
	float lower;
	float upper;

	tkd::CmdParser parser( argv[0], "Calculate mean between bounds; ( included voxels: <= lower, >= upper ) [supports 2D,3D and 4D]" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "mask" ) -> AddAlias( "m" ) -> SetInput( "<string>" ) -> SetDescription( "Mask File" ) -> SetRequired( true ) -> SetMinMax(
			1, 1 );

	parser.AddArgument( lower, "lower" ) -> AddAlias( "l" ) -> SetInput( "<float>" ) -> SetDescription( "Lower bound" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( upper, "upper" ) -> AddAlias( "u" ) -> SetInput( "<float>" ) -> SetDescription( "Upper bound" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Mean mean = Mean();

	mean.run( input, output, lower, upper );
}

