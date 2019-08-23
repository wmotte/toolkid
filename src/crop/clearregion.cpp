/*
 * clearregion.cpp
 *
 *  Created on: Jun 26, 2009
 *      Author: wim
 */
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

class ClearRegion {

public:
	void run( const std::string& input, const std::string& output, const std::vector< int >& index, const std::vector< int >& size ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 3 ) {
			process( input, output, index, size );
		} else {
			std::cerr << "Number of dimensions 3!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:
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

	void process( const std::string& input, const std::string& output, const std::vector< int >& index, const std::vector< int >& size ) {

		// sanity check...
		if ( ( index.size() != 3 ) && ( size.size() != 3 ) ) {
			std::cerr << "Index and size should be 3!" << std::endl;
			exit( EXIT_FAILURE );
		}

		typedef itk::Image< float, 3 > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

		// get input image...
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		ImageType::Pointer image = reader -> GetOutput();
		reader -> Update();

		ImageType::RegionType clearRegion = image -> GetLargestPossibleRegion();
		ImageType::IndexType clearIndex = clearRegion.GetIndex();
		ImageType::SizeType clearSize = clearRegion.GetSize();

		// TODO

		for ( unsigned int i = 0; i < size.size(); i++ ) {

			if ( ( (unsigned) ( size[i] + index[i] ) ) <= clearSize[i] ) {
				clearSize[i] = size[i];
				clearIndex[i] = index[i];
			} else {
				std::cerr << "Clear region out of bounds!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		clearRegion.SetSize( clearSize );
		clearRegion.SetIndex( clearIndex );

		// find index and size to crop...
		IteratorType it( image, clearRegion );

		// set all voxels in clearRegion
		for ( ; !it.IsAtEnd(); ++it ) {
			image -> SetPixel( it.GetIndex(), 0 );
		}

		// write output...
		WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( output );
		writer -> SetInput( image );

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
	std::vector< int > index;
	std::vector< int > size;

	tkd::CmdParser parser( argv[0], "Clear region in image. (supports 3D only)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( index, "index" ) -> AddAlias( "r" ) -> SetInput( "<ints>" ) -> SetDescription( "Index of clearance area (x,y,z)" ) -> SetRequired(
			true ) -> SetMinMax( 3, 3 );

	parser.AddArgument( size, "size" ) -> AddAlias( "s" ) -> SetInput( "<ints>" ) -> SetDescription( "Size of clearance area (x,y,z)" ) -> SetRequired(
			true ) -> SetMinMax( 3, 3 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ClearRegion clearRegion = ClearRegion();

	clearRegion.run( input, output, index, size );
}

