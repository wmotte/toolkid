/*
 * cog.cpp
 *
 *  Created on: Jun 10, 2009
 *      Author: wim
 */

#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

/**
 * In this class the center of gravity is calcaulated.
 *
 * Definition taken from wikipedia (http://en.wikipedia.org/wiki/Center_of_mass):
 * --------------------------------
 * The center of mass {R} of a system of particles is defined as
 * the average of their positions, {r_i}, weighted by their masses, {m_i}.
 * {R} = sum {m_i} * {r_i} / sum {m_i}.
 * </math>
 */
class COG {

public:

	/**
	 * Write output after resetting origins.
	 */
	void cog( const std::string& inputImageName, bool outputMm ) {

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
			process< pixel, dimension >( inputImageName, outputMm ); \
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
	void process( const std::string& inputImageName, bool outputMm ) {

		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

		// get input image...
		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputImageName );
		typename ImageType::Pointer input = reader -> GetOutput();
		reader -> Update();

		std::vector< TPixel > positions( VDimension, 0 );
		TPixel total_masses = 0;

		typename ImageType::RegionType region = input -> GetLargestPossibleRegion();
		IteratorType it( input, region );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) {

			typename ImageType::IndexType index = it.GetIndex();
			TPixel mass = input -> GetPixel( index );

			for ( unsigned int i = 0; i < VDimension; i++ ) {
				positions[i] += mass * index[i];
			}

			total_masses += mass;
		}

		for ( unsigned int i = 0; i < VDimension; i++ ) {
			positions[i] /= total_masses;
		}

		// convert to MM
		if ( outputMm ) {

			typename ImageType::SpacingType spacing = input -> GetSpacing();

			for ( unsigned int i = 0; i < VDimension; i++ ) {
				positions[i] *= spacing[i];
			}
		}

		// output
		for ( unsigned int i = 0; i < VDimension; i++ ) {
			if ( i != VDimension - 1 ) {
				std::cout << positions[i] << " ";
			} else {
				std::cout << positions[i] << std::endl;
			}
		}
	}
};

/**
 * Set voxel spacing with given values. Update Input.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string inputImageName;
	bool outputMm = false;

	tkd::CmdParser parser( argv[0], "Get Center of Gravity (supports 2D, 3D and 4D images. "
			"Origin is not taken into account!)" );

	parser.AddArgument( inputImageName, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( outputMm, "millimeters" ) -> AddAlias( "mm" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Output center of gravity in mm instead of voxels" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	COG cog = COG();

	cog.cog( inputImageName, outputMm );
}

