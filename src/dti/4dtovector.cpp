#include "tkdCmdParser.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkVector.h"
#include "itkImageFileWriter.h"

int main( int argc, char ** argv )
{
	tkd::CmdParser p( "4dtovector", "Convert 4D image to 3D vector image" );

	std::string inputFileName;
	std::string outputFileName;
	int peak = 0;
	bool keepNaN = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input 4D image" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output filename" ) ->SetRequired( true );

	p.AddArgument( peak, "peak" ) ->AddAlias( "p" ) ->SetDescription( "Peak number to extract (default: 0)" );

	p.AddArgument( keepNaN, "keep-nan" ) ->AddAlias( "kn" ) ->SetDescription( "Do not zero NaN values" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	typedef float PixelType;
	typedef itk::Vector< PixelType, 3 > VectorPixelType;
	typedef itk::Image< PixelType, 4 > InputImageType;
	typedef itk::Image< VectorPixelType, 3 > OutputImageType;
	typedef itk::ImageFileReader< InputImageType > ReaderType;
	typedef itk::ImageFileWriter< OutputImageType > WriterType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inputFileName.c_str() );
	reader->Update();

	InputImageType::Pointer input = reader->GetOutput();

	InputImageType::RegionType inputRegion = input->GetLargestPossibleRegion();
	InputImageType::SizeType inputSize = inputRegion.GetSize();
	InputImageType::IndexType inputIndex = inputRegion.GetIndex();
	InputImageType::SpacingType inputSpacing = input->GetSpacing();
	InputImageType::PointType inputOrigin = input->GetOrigin();
	InputImageType::DirectionType inputDirection = input->GetDirection();

	OutputImageType::RegionType outputRegion;
	OutputImageType::SizeType outputSize;
	OutputImageType::SpacingType outputSpacing;
	OutputImageType::PointType outputOrigin;
	OutputImageType::DirectionType outputDirection;

	for ( int i = 0; i < 3; ++i )
	{
		outputSize[i] = inputSize[i];
		outputOrigin[i] = inputOrigin[i];
		outputSpacing[i] = inputSpacing[i];

		for ( int j = 0; j < 3; ++j )
		{
			outputDirection( i, j ) = inputDirection( i, j );
		}
	}

	outputRegion.SetSize( outputSize );

	OutputImageType::Pointer output = OutputImageType::New();
	output->SetRegions( outputRegion );
	output->SetSpacing( outputSpacing );
	output->SetOrigin( outputOrigin );
	output->SetDirection( outputDirection );
	output->Allocate();
	output->FillBuffer( itk::NumericTraits< VectorPixelType >::Zero );

	typedef itk::ImageRegionConstIterator< InputImageType > InputIterator;
	typedef itk::ImageRegionIterator< OutputImageType > OutputIterator;

	std::vector< InputIterator > inputIt;
	for ( int i = 0; i < 3; ++i )
	{
		inputIndex[3] = peak * 3 + i;
		inputSize[3] = 1;
		inputRegion.SetSize( inputSize );
		inputRegion.SetIndex( inputIndex );
		inputIt.push_back( InputIterator( input, inputRegion ) );
		inputIt[i].GoToBegin();
	}

	for ( OutputIterator outputIt( output, outputRegion ); !outputIt.IsAtEnd(); ++outputIt )
	{
		VectorPixelType vector;

		bool nan = false;

		for ( int i = 0; i < 3; ++i )
		{
			InputIterator& it = inputIt[i];
			nan = nan || ( it.Value() != it.Value() );
			vector[i] = it.Value();
			++it;
		}

		if ( !nan || keepNaN )
		{
			outputIt.Set( vector );
		}
	}

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outputFileName.c_str() );
	writer->SetInput( output );
	writer->Update();

	return 0;
}
