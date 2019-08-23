#include "dtifit.h"
#include "tkdCmdParser.h"
#include "procparser.h"
#include "vnl/vnl_matrix.h"
#include "itkImageFileReader.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "dtifit", "Process DTI data" );

  std::string inputFileName;
  std::string fidFileName;
  std::string outputFileName;
  std::string outputExtension = "nii.gz";
  bool mirrorGradientTable = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input 4D image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output filename base" )
    ->SetRequired( true );

  p.AddArgument( outputExtension, "extension" )
    ->AddAlias( "e" )
    ->SetDescription( "Output filename extension (default: nii.gz)" );

  p.AddArgument( fidFileName, "fid" )
    ->AddAlias( "f" )
    ->SetDescription( "FID filename" )
    ->SetRequired( true );

  p.AddArgument( mirrorGradientTable, "mirror" )
    ->AddAlias( "m" )
    ->SetDescription( "Add negative directions to gradient table" )
    ->SetRequired( false );

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

  typedef float PixelType;
  typedef itk::Image< PixelType, 4 > SeriesType;
  typedef itk::ImageFileReader< SeriesType > ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  dtifit::DTIFit::Pointer fit = dtifit::DTIFit::New();
  fit->Run( reader->GetOutput(), pp, outputFileName, outputExtension, mirrorGradientTable );

  return 0;
}
