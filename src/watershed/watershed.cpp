/*
 * watershed.cpp
 *
 *  Created on: Jun 23, 2009
 *      Author: wim
 */
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorCastImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"

#include "tkdCmdParser.h"

class Watershed {

public:
	/**
	 * Run segementation.
	 */
	void run( const std::string& input, const std::string& output, double level, double threshold, bool pca, unsigned int iterations, double conductance, double timeStep ) {

		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 2 ) {
			segment( input, output, level, threshold, pca, iterations, conductance, timeStep );
		} else {
			std::cerr << "Number of dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	/**
	 * Segment.
	 */
	void segment( const std::string& input, const std::string& output, double level, double threshold,
			bool pca, unsigned int iterations, double conductance, double timeStep ) {

		typedef itk::RGBPixel< unsigned char > RGBPixelType;
		typedef itk::Image< RGBPixelType, 2 > RGBImageType;
		typedef itk::Vector< float, 3 > VectorPixelType;
		typedef itk::Image< VectorPixelType, 2 > VectorImageType;
		typedef itk::Image< unsigned long, 2 > LabeledImageType;
		typedef itk::Image< float, 2 > ScalarImageType;
		typedef itk::ImageFileReader< RGBImageType > FileReaderType;
		typedef itk::VectorCastImageFilter< RGBImageType, VectorImageType > CastFilterType;
		typedef itk::VectorGradientAnisotropicDiffusionImageFilter< VectorImageType, VectorImageType > DiffusionFilterType;
		typedef itk::VectorGradientMagnitudeImageFilter< VectorImageType > GradientMagnitudeFilterType;
		typedef itk::WatershedImageFilter< ScalarImageType > WatershedFilterType;
		typedef itk::ImageFileWriter< RGBImageType > FileWriterType;
		typedef itk::Functor::ScalarToRGBPixelFunctor< unsigned long > ColorMapFunctorType;
		typedef itk::UnaryFunctorImageFilter< LabeledImageType, RGBImageType, ColorMapFunctorType > ColorMapFilterType;
		ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();

		FileReaderType::Pointer reader = FileReaderType::New();
		reader -> SetFileName( input );

		CastFilterType::Pointer caster = CastFilterType::New();

		DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
		diffusion->SetNumberOfIterations( iterations );
		diffusion->SetConductanceParameter( conductance );
		diffusion->SetTimeStep( timeStep ); // timestep...

		GradientMagnitudeFilterType::Pointer gradient = GradientMagnitudeFilterType::New();

		if ( pca ) {
			gradient -> SetUsePrincipleComponentsOn();
		} else {
			gradient -> SetUsePrincipleComponentsOff();
		}

		// There are two parameters that control watershed depth, the lower thresholding of the input.
		// Both parameters are set as a percentage (0.0 - 1.0) of the maximum depth in the input image.
		if ( ( level < 0.0 ) || ( level > 1.0 ) ) {
			std::cerr << "Level parameter, (" << level << ") should be between 0.0 and 1.0!" << std::endl;
			exit( EXIT_FAILURE );
		} else if ( ( threshold < 0.0 ) || ( threshold > 1.0 ) ) {
			std::cerr << "Level parameter, (" << threshold << ") should be between 0.0 and 1.0!" << std::endl;
			exit( EXIT_FAILURE );
		}

		WatershedFilterType::Pointer watershed = WatershedFilterType::New();
		watershed -> SetLevel( level );
		watershed -> SetThreshold( threshold );

		FileWriterType::Pointer writer = FileWriterType::New();
		writer -> SetFileName( output );

		caster->SetInput( reader->GetOutput() );
		diffusion->SetInput( caster->GetOutput() );
		gradient->SetInput( diffusion->GetOutput() );
		watershed->SetInput( gradient->GetOutput() );
		colormapper->SetInput( watershed->GetOutput() );
		writer->SetInput( colormapper->GetOutput() );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
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

};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;
	std::string output;
	double level = 0.0;
	double threshold = 0.0;
	bool pca = true;
	int iterations = 1;
	double conductance = 0.0;
	double timeStep = 0.125;

	tkd::CmdParser parser( argv[0], "Watershed Atlas Image. (supports 2D only)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Image File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( level, "level" ) -> AddAlias( "l" ) -> SetInput( "<double>" ) -> SetDescription(
			"Watershed depth (flood level: default 0.0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( threshold, "threshold" ) -> AddAlias( "t" ) -> SetInput( "<double>" ) -> SetDescription(
			"Lower threshold (default 0.0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( pca, "pca" ) -> AddAlias( "p" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Use Principal Component Analysis to construct Gradient Magnitude (default true)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( iterations, "iterations" ) -> AddAlias( "it" ) -> SetInput( "<uint>" ) -> SetDescription(
				"Vector Gradient Anisotropic Diffusion Iterations (smoothing) (default 1)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( conductance, "conductance" ) -> AddAlias( "c" ) -> SetInput( "<double>" ) -> SetDescription(
					"Vector Gradient Anisotropic Diffusion Conductance (smoothing) (default 0.0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( timeStep, "time-step" ) -> AddAlias( "s" ) -> SetInput( "<double>" ) -> SetDescription(
						"Vector Gradient Anisotropic Diffusion Time Step (smoothing) (default 0.125)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );


	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Watershed watershed = Watershed();

	watershed.run( input, output, level, threshold, pca, (unsigned) iterations, conductance, timeStep );
}
