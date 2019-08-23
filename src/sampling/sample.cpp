#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkIdentityTransform.h"
#include "tkdCmdParser.h"
#include <string>
#include <itkNearestNeighborInterpolateImageFunction.h>

class Sample
{

public:
	void run( const std::string& input, const std::string& reference, const std::string& output, bool mask, bool series )
	{

		if ( !series )
		{
			typedef itk::Image< float, 3 > ImageType;
			typedef itk::ImageFileReader< ImageType > ReaderType;
			typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleType;
			typedef itk::ImageFileWriter< ImageType > WriterType;
			typedef itk::IdentityTransform< double, 3 > TransformType;

			ReaderType::Pointer inputReader = ReaderType::New();
			ReaderType::Pointer referenceReader = ReaderType::New();
			ResampleType::Pointer resample = ResampleType::New();
			WriterType::Pointer writer = WriterType::New();
			TransformType::Pointer transform = TransformType::New();

			inputReader -> SetFileName( input );
			referenceReader -> SetFileName( reference );
			writer -> SetFileName( output );

			inputReader -> Update();
			referenceReader -> Update();

			resample -> SetTransform( transform );
			resample -> SetInput( inputReader -> GetOutput() );
			resample -> SetReferenceImage( referenceReader -> GetOutput() );
			resample -> UseReferenceImageOn();
			resample -> Update();

			if ( mask )
			{
				typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;

				InterpolatorType::Pointer interpolator = InterpolatorType::New();
				resample -> SetInterpolator( interpolator );
			}

			writer -> SetInput( resample -> GetOutput() );

			try
			{
				writer -> Update();
			} catch ( itk::ExceptionObject & exp )
			{
				std::cerr << "ERROR - could not write: " << output << std::endl;
				exit( EXIT_FAILURE );
			}
		}
		else // 4D
		{
			typedef itk::Image< float, 4 > ImageType;
			typedef itk::ImageFileReader< ImageType > ReaderType;
			typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleType;
			typedef itk::ImageFileWriter< ImageType > WriterType;
			typedef itk::IdentityTransform< double, 4 > TransformType;

			ReaderType::Pointer inputReader = ReaderType::New();
			ReaderType::Pointer referenceReader = ReaderType::New();
			ResampleType::Pointer resample = ResampleType::New();
			WriterType::Pointer writer = WriterType::New();
			TransformType::Pointer transform = TransformType::New();

			inputReader -> SetFileName( input );
			referenceReader -> SetFileName( reference );
			writer -> SetFileName( output );

			inputReader -> Update();
			referenceReader -> Update();

			resample -> SetTransform( transform );
			resample -> SetInput( inputReader -> GetOutput() );
			resample -> SetReferenceImage( referenceReader -> GetOutput() );
			resample -> UseReferenceImageOn();
			resample -> Update();

			if ( mask )
			{
				typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;

				InterpolatorType::Pointer interpolator = InterpolatorType::New();
				resample -> SetInterpolator( interpolator );
			}

			writer -> SetInput( resample -> GetOutput() );

			try
			{
				writer -> Update();
			} catch ( itk::ExceptionObject & exp )
			{
				std::cerr << "ERROR - could not write: " << output << std::endl;
				exit( EXIT_FAILURE );
			}
		}
	}
};

/**
 * Sample.
 */
int main( int argc, char ** argv )
{

	// arguments...
	std::string input;
	std::string reference;
	std::string output;
	bool mask = false;
	bool series = false;

	tkd::CmdParser parser( argv[0], "Sample one image to spacing of reference image" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input image" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( reference, "reference" ) -> AddAlias( "r" ) -> SetInput( "<string>" ) -> SetDescription( "Reference image" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output image" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( mask, "mask" ) -> AddAlias( "m" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Nearest neighbor interpolation (default: false)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( series, "series" ) -> AddAlias( "s" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Series input (default: false)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Sample sample = Sample();
	sample.run( input, reference, output, mask, series );
}

