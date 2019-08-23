/*
 * histogram.cpp
 *
 *  Created on: Sep 28, 2009
 *      Author: wim
 */
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkScalarImageToListAdaptor.h"
#include "itkListSampleToHistogramGenerator.h"

/**
 * Output cumulative image histogram.
 */
class CumulativeHistogram {

public:

	void run( const std::string& input, const std::string& output, unsigned int bins, float minValue, float maxValue, bool percentage ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 2 ) {
			process2D( input, output, bins, minValue, maxValue, percentage );
		} else if ( dims == 3 ) {
			process3D( input, output, bins, minValue, maxValue, percentage );
		} else {
			std::cerr << "Number of dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

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

	/**
	 * Process 2D input image.
	 */
	void process2D( const std::string& input, const std::string& output, unsigned int bins, float minValue, float maxValue, bool percentage ) {

		typedef float PixelType;
		typedef itk::Image< PixelType, 2 > Image2DType;
		typedef itk::ImageFileReader< Image2DType > Reader2DType;

		typedef itk::ImageFileWriter< Image2DType > Writer2DType;
		typedef itk::Statistics::ScalarImageToListAdaptor< Image2DType > Adaptor2DType;
		typedef PixelType HistogramMeasurementType;
		typedef itk::Statistics::ListSampleToHistogramGenerator< Adaptor2DType, HistogramMeasurementType > GeneratorType;
		typedef GeneratorType::HistogramType HistogramType;

		// read...
		Reader2DType::Pointer reader = Reader2DType::New();
		reader -> SetFileName( input );
		reader -> Update();

		// adaptor...
		Adaptor2DType::Pointer adaptor = Adaptor2DType::New();
		adaptor->SetImage( reader->GetOutput() );

		// generator...
		GeneratorType::Pointer generator = GeneratorType::New();
		HistogramType::SizeType size;

		size.Fill( bins );

		generator->SetListSample( adaptor );
		generator->SetNumberOfBins( size );
		generator->SetMarginalScale( 10.0 );

		HistogramType::MeasurementVectorType min;
		HistogramType::MeasurementVectorType max;

		min.Fill( minValue );
		max.Fill( maxValue );

		generator->SetHistogramMin( min );
		generator->SetHistogramMax( max );

		generator->Update();

		HistogramType::ConstPointer histogram = generator->GetOutput();

		const unsigned int histogramSize = histogram->Size();

		// write output...
		std::ofstream filestream;
		filestream.open( output.c_str() );
		if ( filestream.fail() ) {
			std::cerr << "Not able to write to: " << output << std::endl;
			exit( EXIT_FAILURE );
		} else {

			// get total frequencies
			double totalFrequencies = 0;
			for ( unsigned int i = 0; i < histogramSize; i++ ) {
				totalFrequencies += histogram->GetFrequency( i, 0 );
			}

			// get cumulative frequency
			for ( unsigned int i = 0; i < histogramSize; i++ ) {
				double cumFrequency = 0;
				for ( unsigned int j = i; j < histogramSize; j++ ) {
					cumFrequency += histogram->GetFrequency( j, 0 );
				}

				if ( percentage ) {
					filestream << ( ( cumFrequency / totalFrequencies ) * 100 ) << std::endl;
				} else {
					filestream << ( cumFrequency / totalFrequencies ) << std::endl;
				}
			}
		}
		filestream.close();
	}

	/**
	 * Process 3D input image.
	 */
	void process3D( const std::string& input, const std::string& output, unsigned int bins, float minValue, float maxValue, bool percentage ) {

		typedef float PixelType;
		typedef itk::Image< PixelType, 3 > Image3DType;
		typedef itk::ImageFileReader< Image3DType > Reader3DType;

		typedef itk::ImageFileWriter< Image3DType > Writer3DType;
		typedef itk::Statistics::ScalarImageToListAdaptor< Image3DType > Adaptor3DType;
		typedef PixelType HistogramMeasurementType;
		typedef itk::Statistics::ListSampleToHistogramGenerator< Adaptor3DType, HistogramMeasurementType > GeneratorType;
		typedef GeneratorType::HistogramType HistogramType;

		// read...
		Reader3DType::Pointer reader = Reader3DType::New();
		reader -> SetFileName( input );
		reader -> Update();

		// adaptor...
		Adaptor3DType::Pointer adaptor = Adaptor3DType::New();
		adaptor->SetImage( reader->GetOutput() );

		// generator...
		GeneratorType::Pointer generator = GeneratorType::New();
		HistogramType::SizeType size;

		size.Fill( bins );

		generator->SetListSample( adaptor );
		generator->SetNumberOfBins( size );
		generator->SetMarginalScale( 10.0 );

		HistogramType::MeasurementVectorType min;
		HistogramType::MeasurementVectorType max;

		min.Fill( minValue );
		max.Fill( maxValue );

		generator->SetHistogramMin( min );
		generator->SetHistogramMax( max );

		generator->Update();

		HistogramType::ConstPointer histogram = generator->GetOutput();

		const unsigned int histogramSize = histogram->Size();

		// write output...
		std::ofstream filestream;
		filestream.open( output.c_str() );
		if ( filestream.fail() ) {
			std::cerr << "Not able to write to: " << output << std::endl;
			exit( EXIT_FAILURE );
		} else {

			// get total frequencies
			double totalFrequencies = 0;
			for ( unsigned int i = 0; i < histogramSize; i++ ) {
				totalFrequencies += histogram->GetFrequency( i, 0 );
			}

			// get cumulative frequency
			for ( unsigned int i = 0; i < histogramSize; i++ ) {
				double cumFrequency = 0;
				for ( unsigned int j = i; j < histogramSize; j++ ) {
					cumFrequency += histogram->GetFrequency( j, 0 );
				}

				if ( percentage ) {
					filestream << ( ( cumFrequency / totalFrequencies ) * 100 ) << std::endl;
				} else {
					filestream << ( cumFrequency / totalFrequencies ) << std::endl;
				}
			}
		}
		filestream.close();
	}

};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;
	std::string output;
	int bins = 1000;
	float minValue = 0;
	float maxValue = 100;
	bool percentage = false;

	tkd::CmdParser parser( argv[0], "Calculate histogram of input image." );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Text File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( bins, "bins" ) -> AddAlias( "b" ) -> SetInput( "<int>" ) -> SetDescription( "Number of bins (default: 1000)" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( minValue, "min" ) -> AddAlias( "l" ) -> SetInput( "<float>" ) -> SetDescription( "Min value (default: 0)" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( maxValue, "max" ) -> AddAlias( "u" ) -> SetInput( "<float>" ) -> SetDescription( "Max value (default: 100)" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( percentage, "percentage" ) -> AddAlias( "p" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Convert absolute frequencies to percentage of total frequencies (default: false)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	CumulativeHistogram histogram = CumulativeHistogram();

	histogram.run( input, output, (unsigned) bins, minValue, maxValue, percentage );
}
