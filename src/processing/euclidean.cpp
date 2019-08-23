
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "tkdCmdParser.h"


/**
 * Euclidean distance to cog label image.
 */
class Euclidean
{

public:

	typedef float PixelType;
	typedef itk::Point< PixelType, 3 > PointType;
	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ConstIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

	/**
	 * run.
	 */
	void run( const std::string& mask, const std::string& label, const std::string& output )
	{
		process( mask, label, output );

		exit( EXIT_SUCCESS );
	}

	/**
	 * Return center-of-gravity for 3D volume.
	 */
	PointType getCog( const ImageType::ConstPointer volume )
	{
		PixelType posx = 0;
		PixelType posy = 0;
		PixelType posz = 0;
		PixelType total_masses = 0;

		ImageType::RegionType region = volume -> GetLargestPossibleRegion();

		ImageType::SpacingType spacing = volume -> GetSpacing();

		ConstIteratorType it( volume, region );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			ImageType::IndexType index = it.GetIndex();
			PixelType mass = volume -> GetPixel( index );

			posx += mass * index[0];
			posy += mass * index[1];
			posz += mass * index[2];

			total_masses += mass;
		}

		PointType point;
		point[0] = ( posx * spacing[0] ) / total_masses;
		point[1] = ( posy * spacing[1] ) / total_masses;
		point[2] = ( posz * spacing[2] ) / total_masses;

		return point;
	}


	/**
	 * Process.
	 */
	void process( const std::string& maskFileName, const std::string& labelFileName, const std::string& outputFileName )
	{
		// read mask...
		ReaderType::Pointer maskReader = ReaderType::New();
		maskReader -> SetFileName( maskFileName );
		maskReader -> Update();

		// read label...
		ReaderType::Pointer labelReader = ReaderType::New();
		labelReader -> SetFileName( labelFileName );
		labelReader -> Update();

		ImageType::ConstPointer mask = maskReader -> GetOutput();
		ImageType::ConstPointer label = labelReader -> GetOutput();

		ImageType::RegionType maskRegion = mask -> GetLargestPossibleRegion();
		ImageType::RegionType labelRegion = label -> GetLargestPossibleRegion();

		ImageType::SizeType maskSize = maskRegion.GetSize();
		ImageType::SizeType labelSize = labelRegion.GetSize();

		ImageType::SpacingType spacing = mask -> GetSpacing();

		// allocate new image
		ImageType::Pointer output = ImageType::New();
		output->CopyInformation( label );
		output->SetRegions( label->GetLargestPossibleRegion() );
		output->Allocate();
		output->FillBuffer( 0 );

		ImageType::RegionType outputRegion = output->GetLargestPossibleRegion();

		// size check...
		for( unsigned int i = 0; i < 3; i++ )
		{
			if( labelSize[i] != maskSize[i] )
			{
				std::cerr << "Label and Mask Size don't match!" << labelSize << " != " << maskSize << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		// *******************************

		PointType cog = getCog( label );

		std::cout << "cog: " << cog << std::endl;

		std::vector< PixelType > includedValues;
		ConstIteratorType mit( mask, maskRegion );
		IteratorType ot( output, outputRegion );

		for ( ; ot.IsAtEnd(), !mit.IsAtEnd(); ++mit, ++ot )
		{
			ImageType::IndexType index = mit.GetIndex();
			PointType point;

			point[ 0 ] = spacing[ 0 ] * index[ 0 ];
			point[ 1 ] = spacing[ 1 ] * index[ 1 ];
			point[ 2 ] = spacing[ 2 ] * index[ 2 ];

			if ( mask -> GetPixel( index ) != 0 )
			{
				ot.Set( point.EuclideanDistanceTo( cog  ) );
			}
		}

		WriterType::Pointer writer = WriterType::New();
		writer->SetInput( output );
		writer->SetFileName( outputFileName );
		writer->Update();



	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::string mask;
	std::string label;
	std::string output;

	tkd::CmdParser parser( argv[0], "Calculate euclidean distance to cog of label image" );

	parser.AddArgument( mask, "mask" ) -> AddAlias( "m" ) -> SetInput( "<string>" ) -> SetDescription( "Mask image" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );
	parser.AddArgument( label, "label" ) -> AddAlias( "l" ) -> SetInput( "<string>" ) -> SetDescription( "Label image" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );
	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output image" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Euclidean r = Euclidean();

	r.run( mask, label, output );
}

