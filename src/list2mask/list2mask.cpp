#include "itkImageRegionConstIterator.h"
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

/**
 * In this class a text-file containing values for each voxel in a mask is plotted back on the mask.
 * This is useful when voxel-wise graph analysis methods are use to calculate matrix properties (e.g. clustering).
 * This class is able to project this matrix properties back to the original voxels, given the original mask.
 */
class List2Mask {



public:
	typedef float PixelType;
	typedef itk::Image< PixelType, 3> ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	static const PixelType BACKGROUND_VALUE = 0;


	// get vector...
	void textToVector( const std::string& inputTextName, std::vector<int>& vector ) {

		int temp;
		std::ifstream istr;

		istr.open( inputTextName.c_str() );
		if ( istr.fail() ) {
			std::cerr << "Failed to read: " << inputTextName << std::endl;
		}

		while ( istr >> temp ) {
			vector.push_back( temp );
		}
	}

	// get input image...
	void getImage( const std::string& inputImageName, ImageType::Pointer& inputImage ) {
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputImageName );
		reader -> Update();
		inputImage = reader -> GetOutput();

		reader -> Delete();
	}

	// get output image...
	void vectorToImage( const std::vector<int>& vector, const ImageType::Pointer& inputImage, ImageType::Pointer& outputImage ) {

		outputImage = ImageType::New();
		outputImage -> SetRegions( inputImage -> GetLargestPossibleRegion() );
		outputImage -> SetOrigin( inputImage -> GetOrigin() );
		outputImage -> SetSpacing( inputImage -> GetSpacing() );
		outputImage -> Allocate();


		itk::ImageRegionConstIterator< ImageType > iterator( inputImage, inputImage -> GetLargestPossibleRegion() );

		unsigned int position = 0;

		for( iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator ) {

			if ( iterator.Get() != BACKGROUND_VALUE ) {

				iterator.GetIndex();
				if ( vector.size() == position ) {
					std::cerr << "ERROR: Vector size ("<< vector.size() << ") is smaller than number of nonzero voxels!" << std::endl;
				}

				outputImage -> SetPixel( iterator.GetIndex(), vector[ position ] );
				position++;
			}
		}
	}

	// write output image...
	void writeImage( const std::string& outputImageName, const ImageType::Pointer& outputImage ) {
		WriterType::Pointer writer = WriterType::New();
		writer -> SetInput( outputImage );
		writer -> SetFileName( outputImageName );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& exception ) {
			std::cerr << "Error writing: " << outputImageName << std::endl;
			std::cerr << exception.GetDescription() << std::endl;
		}
	}

};

int main( int argc, char * argv[] ) {

	// arguments...
	std::string inputTextName;
	std::string inputImageName;
	std::string outputImageName;

	tkd::CmdParser parser( argv[0], "Map list (txt) to Image (3D)" );

	parser.AddArgument( inputTextName, "text" )
		-> AddAlias( "t" )
		-> SetInput( "<string>")
		-> SetDescription( "Input Text File" )
		-> SetRequired( true )
		-> SetMinMax( 1, 1 );

	parser.AddArgument( inputImageName, "input" )
		-> AddAlias( "i" )
		-> SetInput( "<string>" )
		-> SetDescription( "Input Image File" )
		-> SetRequired( true )
		-> SetMinMax( 1, 1 );

	parser.AddArgument( outputImageName, "output")
		-> AddAlias( "o" )
		-> SetInput( "<string>" )
		-> SetDescription( "Output Image File" )
		-> SetRequired( true )
		-> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
	    return EXIT_FAILURE;
	}


	List2Mask list2Mask;

	std::vector<int> vector;
	List2Mask::ImageType::Pointer inputImage;
	List2Mask::ImageType::Pointer outputImage;

	// get vector...
	list2Mask.textToVector( inputTextName, vector );

	// get input image...
	list2Mask.getImage( inputImageName, inputImage );

	// get output image...
	list2Mask.vectorToImage( vector, inputImage, outputImage );

	// write output image...
	list2Mask.writeImage( outputImageName, outputImage );

}





