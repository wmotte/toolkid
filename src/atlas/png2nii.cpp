/*
 * watershed.cpp
 *
 *  Created on: Jun 29, 2009
 *      Author: wim
 */
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNormalizeImageFilter.h"

/**
 * Convert multiple pngs to 3d volume.
 */

class Png2Nii {

public:
	void run( const std::vector< std::string >& inputs, const std::string& output ) {

		// get dimensions...
		if ( dimensionsOke( inputs ) ) {
			process( inputs, output );
		} else {
			std::cerr << "Input images have different dimensions!" << std::endl;
			exit( EXIT_FAILURE );
		}

	}

protected:

	/**
	 * Process 2D inputs...
	 */
	void process( const std::vector< std::string >& inputs, const std::string& outputFileName ) {
		// sanity check...
		if ( inputs.empty() ) {
			std::cerr << "No inputs" << std::endl;
			exit( EXIT_FAILURE );
		}

		typedef float PixelType;
		typedef itk::Image< PixelType, 3 > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ConstIteratorType;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
		typedef itk::NormalizeImageFilter<ImageType, ImageType> NormalizeImageFilterType;

		// A. First get size from intput files...
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputs[0] );
		ImageType::ConstPointer slice = reader -> GetOutput();
		reader -> Update();
		ImageType::RegionType region = slice -> GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		ImageType::IndexType index = region.GetIndex();
		size[2] = inputs.size();
		region.SetSize( size );

		// B. Next, allocate output...
		ImageType::Pointer output = ImageType::New();
		output -> SetRegions( region );
		output -> Allocate();

		NormalizeImageFilterType::Pointer filter = NormalizeImageFilterType::New();

		// C. Read each slice and insert in output volume...
		for ( unsigned int i = 0; i < inputs.size(); i++ ) {
			reader -> SetFileName( inputs[i] );
			filter -> SetInput( reader -> GetOutput() );
			slice = filter -> GetOutput();
			filter -> Update();
			ImageType::RegionType sliceRegion = slice -> GetLargestPossibleRegion();
			ConstIteratorType sit( slice, sliceRegion );

			// for each slice pixel...
			for ( sit.GoToBegin(); !sit.IsAtEnd(); ++sit ) {
				ImageType::IndexType index = sit.GetIndex();
				PixelType pixel = slice -> GetPixel( index );

				// update index to 3D represenation...
				index[2] = i;
				// ...and inset in 3d volume.
				output -> SetPixel( index, pixel );
			}
		}

		WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( outputFileName );
		writer -> SetInput( output );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing to: " << outputFileName << std::endl;
			std::cerr << e.GetDescription() << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Return false if one of the inputs has dims != 2.
	 */
	bool dimensionsOke( const std::vector< std::string > inputs ) {
		// check every input file...
		for ( unsigned int i = 0; i < inputs.size(); i++ ) {
			std::string inputFileName = inputs[i];
			itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
					itk::ImageIOFactory::ReadMode );
			if ( !io ) {
				std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
				exit( EXIT_FAILURE );
			} else {
				io -> SetFileName( inputFileName );
				io -> ReadImageInformation();
				unsigned int dims = io -> GetNumberOfDimensions();
				if ( dims != 2 ) {
					io = 0;
					return false;
				}
			}
		}
		return true;
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::vector< std::string > inputs;
	std::string output;

	tkd::CmdParser parser( argv[0], "Convert 2D Images to 3D ITK Image." );

	parser.AddArgument( inputs, "inputs" ) -> AddAlias( "i" ) -> SetInput( "<strings>" ) -> SetDescription( "Input 2D Files (e.g. pngs)" ) -> SetRequired(
			true ) -> SetMinMax( 1, 100000 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output 3D File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Png2Nii png2nii = Png2Nii();
	png2nii.run( inputs, output );
}

