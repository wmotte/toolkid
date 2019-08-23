#include "dtifit.h"
#include "tkdCmdParser.h"
#include "procparser.h"
#include "vnl/vnl_matrix.h"
#include "itkImageFileReader.h"

int main( int argc, char ** argv )
{
	tkd::CmdParser p( "dti", "Process DTI data (linear or nonlinear)" );

	std::string inputFileName;
	std::string maskFileName;
	std::string fidFileName;
	std::string outputFileName;
	bool mirrorGradientTable = false;
	int algorithm = 3;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input 4D image" ) ->SetRequired( true );

	p.AddArgument( maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask" );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output filename base" ) ->SetRequired( true );

	p.AddArgument( fidFileName, "fid" ) ->AddAlias( "f" ) ->SetDescription( "FID filename" ) ->SetRequired( true );

	p.AddArgument( mirrorGradientTable, "mirror" ) ->AddAlias( "m" ) ->SetDescription( "Add negative directions to gradient table" ) ->SetRequired( false );

	p.AddArgument( algorithm, "algorithm" ) ->AddAlias( "a" ) ->SetDescription( "Fitting algorithm. 1 => Weighted Linear, 2 => Constrained Weighted Linear, 3 => Constrained Weighted NonLinear (default = 3)" ) ->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	Procparser pp;
	if ( !pp.Parse( fidFileName ) )
	{
		std::cout << "Could not read procpar from " << fidFileName << std::endl;
		return -1;
	}

	typedef double PixelType;
	typedef itk::Image< PixelType, 4 > SeriesType;
	typedef itk::ImageFileReader< SeriesType > ReaderType;
	typedef itk::Image< PixelType, 3 > MaskType;
	typedef itk::ImageFileReader< MaskType > MaskReaderType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inputFileName.c_str() );
	reader->Update();

	dtifit::DTIFit::Pointer fit = dtifit::DTIFit::New();
	fit->Run( reader->GetOutput(), maskFileName, pp, outputFileName, mirrorGradientTable, algorithm );

	return EXIT_SUCCESS;
}
