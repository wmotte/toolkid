#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"

#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "tkdCmdParser.h"

/**
 * Date: 16-11-2009
 */
class CreateRegular {

public:

	typedef double PixelType;
	typedef itk::Image< PixelType, 2 > ImageType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef itk::ImageFileWriter< ImageType > WriterType;

	/**
	 * Run.
	 */
	void run( const std::string& output, unsigned int size, unsigned int neighbors )
	{
		WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( output.c_str() );
		writer->SetInput( GetMatrix( size, neighbors ) );

		try {
			writer->Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}

	}

protected:

	/**
	 * Return regular matrix.
	 */
	ImageType::Pointer GetMatrix( unsigned int matrixSize, unsigned int neighbors )
	{
		ImageType::Pointer outputImage = ImageType::New();
		ImageType::SizeType size;
		ImageType::RegionType region;
		size[ 0 ] = matrixSize;
		size[ 1 ] = matrixSize;

		region.SetSize( size );
		outputImage->SetRegions( region );
		outputImage->Allocate();

		for ( unsigned int i = 0; i < matrixSize; ++i )
		{
			for ( unsigned int j = 0; j < matrixSize; ++j )
			{
				unsigned int distance = abs( i - j );

				// wrap around to get circular network ...
				if ( ( distance <= neighbors && distance > 0 )
					|| (( distance + neighbors ) >= matrixSize ) )
				{
					ImageType::IndexType index;
					index[0] = i;
					index[1] = j;

					outputImage->SetPixel( index, static_cast< PixelType >( 1.0 ) );
				}
			}
		}

		return outputImage;
	}
};

/**
 * Create regular matrix.
 */
int main( int argc, char ** argv ) {
	tkd::CmdParser p( "create_regular", "Construct regular graph." );

	std::string output;
	int size = 50;
	int neighbors = 2;

	p.AddArgument( output, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output matrix file name" ) ->SetRequired(
			true ) -> SetMinMax( 1, 1 );

	p.AddArgument( size, "size" ) ->AddAlias( "s" ) ->SetInput( "uint" ) ->SetDescription( "Number of nodes [ default: 50 ]" );

	p.AddArgument( neighbors, "neighbors" ) ->AddAlias( "n" ) ->SetInput( "uint" ) ->SetDescription( "Number of neighboring nodes [ default: 3 ]" );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	CreateRegular createRegular;

	createRegular.run( output, size, neighbors );

	return EXIT_SUCCESS;
}


