
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"

#include <iomanip>
#include <vector>

/**
 * Max dim.
 */
class MaxDim
{

public:

	/**
	 * Run.
	 */
	void Run( const std::string& input, const std::string& mask, bool mmSize, const std::string& output )
	{

		// get dimensions...
		unsigned int dims = GetDimensions( input );

		if ( dims != 4 )
		{
			std::cerr << "Number of dimensions should be 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		Process( input, mask, mmSize, output );

		exit( EXIT_SUCCESS );
	}

protected:

	typedef float PixelType;
	typedef itk::Image< PixelType, 4 > Image4DType;
	typedef itk::Image< PixelType, 3 > Image3DType;
	typedef itk::ImageFileReader< Image4DType > ReaderType;
	typedef itk::Image< unsigned char, 3 > MaskType;
	typedef itk::ImageFileReader< MaskType > MaskReaderType;
	typedef itk::ImageFileWriter< Image3DType > Image3DWriterType;
	typedef itk::ImageRegionIteratorWithIndex< Image3DType > IteratorType;
	typedef itk::ImageRegionIteratorWithIndex< MaskType > MaskIteratorType;
	typedef itk::ExtractImageFilter< Image4DType, Image3DType > ExtractImageFilterType;

	/**
	 * Check dimensions of inputfile. In case of error, Application is terminated.
	 */
	unsigned int GetDimensions( const std::string& inputFileName )
	{
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io )
		{
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else
		{
			io->SetFileName( inputFileName );
			io->ReadImageInformation();
			return io->GetNumberOfDimensions();
		}
	}

	/**
	 * Process input image.
	 */
	void Process( const std::string& input, const std::string& mask, bool mmSize, const std::string& output )
	{
		// read 4D input
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( input );
		reader->Update();

		// read mask (3D).
		MaskReaderType::Pointer maskReader = MaskReaderType::New();
		maskReader->SetFileName( mask );
		maskReader->Update();

		// for each volume...
		std::vector< MaskType::IndexType > maxPoints;
		std::vector< PixelType > maxValues;

		Image4DType::SizeType size4D = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

		for( unsigned int i = 0; i < size4D[3]; i++ )
		{
			ExtractImageFilterType::Pointer extract = ExtractImageFilterType::New();

			extract->SetInput( reader->GetOutput() );

			Image4DType::RegionType region4D = reader->GetOutput()->GetLargestPossibleRegion();
			Image4DType::IndexType index4D = region4D.GetIndex();
			Image4DType::SizeType size4D = region4D.GetSize();
			size4D[3] = 0;
			index4D[3] = i;
			region4D.SetSize( size4D );
			region4D.SetIndex( index4D );

			extract->SetExtractionRegion( region4D );
			extract->Update();

			MaskType::Pointer maskObject = maskReader->GetOutput();
			Image3DType::Pointer imageObject = extract->GetOutput();

			IteratorType it( extract->GetOutput(), extract->GetOutput()->GetLargestPossibleRegion( ) );
			MaskIteratorType mit( maskReader->GetOutput(), maskReader->GetOutput()->GetLargestPossibleRegion() );

			PixelType maxValue = itk::NumericTraits< PixelType >::min();
			MaskType::IndexType maxPoint;

			for ( it.GoToBegin(), mit.GoToBegin(); !mit.IsAtEnd(); ++mit )
			{
				if ( maskObject->GetPixel( mit.GetIndex() ) != 0 )
				{
					PixelType value = imageObject->GetPixel( mit.GetIndex() );
					if( value > maxValue )
					{
						maxValue = value;
						maxPoint = mit.GetIndex();
					}
				}
			}

			maxPoints.push_back( maxPoint );
			maxValues.push_back( maxValue );
		}

		// print coordinates
		if( mmSize )
			std::cout << "x,y,z,int coordinates of voxel with maximum intensity (mm)." << std::endl;
		else
			std::cout << "x,y,z,int coordinates of voxel with maximum intensity (voxels)." << std::endl;

		for( unsigned int i = 0; i < maxPoints.size(); i++ )
		{
			MaskType::SpacingType spacing = maskReader->GetOutput()->GetSpacing();
			MaskType::PointType origin    = maskReader->GetOutput()->GetOrigin();

			PixelType x = maxPoints.at( i )[0];
			PixelType y = maxPoints.at( i )[1];
			PixelType z = maxPoints.at( i )[2];

			if( mmSize )
			{
				x = ( origin[0] + maxPoints.at( i )[0] * spacing[0] );
				y = ( origin[1] + maxPoints.at( i )[1] * spacing[1] );
				z = ( origin[2] + maxPoints.at( i )[2] * spacing[2] );
			}

			std::cout << std::setprecision( 3 ) << x << ", " << y << ", " << z << ": " <<
			maxValues.at( i ) << std::endl;
		}

		// write max voxels to output...
		if( !output.empty() )
		{
			Image3DType::Pointer newImage = Image3DType::New();
			newImage->CopyInformation( maskReader->GetOutput() );
		    newImage->SetBufferedRegion( maskReader->GetOutput()->GetBufferedRegion() );
		    newImage->SetRequestedRegion( maskReader->GetOutput()->GetRequestedRegion() );

			newImage->Allocate();
			newImage->FillBuffer( 0 );

			// fill 3D volume with max-voxel positions.
			for( unsigned int i = 0; i < maxPoints.size(); i++ )
			{
				newImage->SetPixel( maxPoints.at( i ), static_cast< PixelType >( i + 1. ) );
			}

			Image3DWriterType::Pointer writer = Image3DWriterType::New();
			writer->SetInput( newImage );
			writer->SetFileName( output );

			try
			{
				writer->Update();
			}
			catch( itk::ExceptionObject& e )
			{
				std::cerr << "*** ERROR ***: could not write: " << output << "!" << std::endl;
				std::cerr << e.GetDescription() << std::endl;
				exit( EXIT_FAILURE );
			}
		}

	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::string input;
	std::string output;
	std::string mask;
	bool mmSize = false;

	tkd::CmdParser parser( argv[0], "Calculate coordinates of voxel with maximum intensity." );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Image Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( mask, "mask" ) -> AddAlias( "m" ) -> SetInput( "<string>" ) -> SetDescription( "Mask Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( mmSize, "mm" ) -> AddAlias( "mm" ) -> SetInput( "<bool>" ) -> SetDescription( "Print in mm (default: false)" ) -> SetRequired(
				false );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	MaxDim maxDim = MaxDim();

	maxDim.Run( input, mask, mmSize, output );
}

