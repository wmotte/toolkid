#include <iostream>
#include "mrfitCommon.h"
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "procparser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "t2map", "T2 mapping" );

  std::string inputFileName;
  std::string outputFileName;
  std::string procparFileName;
  std::string outputFileNameA;
  std::string outputFileNameC;
  std::string outputFileNameR2;

  int algorithm = 0;
  double maxT2 = 10.0;
  std::vector< double > echos;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input 4D image (3D x number of echos)" )
    ->SetRequired( true );

  p.AddArgument( procparFileName, "procpar" )
    ->AddAlias( "p" )
    ->SetDescription( "Path to FID/procpar" );

  p.AddArgument( echos, "echos" )
    ->AddAlias( "e" )
    ->SetDescription( "Echo times (in seconds)" );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output T2 map" )
    ->SetRequired( true );

  p.AddArgument( outputFileNameA, "output-A" )
    ->AddAlias( "oA" )
    ->SetDescription( "Output map of fitted constant A" );

  p.AddArgument( outputFileNameC, "output-C" )
    ->AddAlias( "oC" )
    ->SetDescription( "Output map of fitted constant C" );

  p.AddArgument( outputFileNameR2, "output-R2" )
    ->AddAlias( "oR2" )
    ->SetDescription( "Output map of R^2 fitting values" );

  p.AddArgument( algorithm, "algorithm" )
    ->AddAlias( "a" )
    ->SetDescription( "Fitting algorithm: 0=LINEAR, 1=NON_LINEAR, 2=NON_LINEAR_WITH_CONSTANT (default: 0)" );

  p.AddArgument( maxT2, "maximum-T2" )
    ->AddAlias( "max" )
    ->SetDescription( "Maximum T2 (default: 10)" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  Procparser pp;

  if ( echos.size() == 0 )
    {
    if ( procparFileName == "" )
      {
      std::cerr << "Either provide the path to the FID's procpar, or supply the echo times" << std::endl;
      return -1;
      }

    pp.Parse( procparFileName );
    if ( pp.Has( "TE" ) )
      {
      for( int i = 0; i < pp.GetSize( "TE" ); ++i )
        {
        echos.push_back( pp.GetAs< double >( "TE", i ) / 1000. );
        }
      }
    else
      {
      for( int i = 1; i <= pp.GetAs< int >( "ne" ); ++i )
        {
        echos.push_back( pp.GetAs< double >( "te" ) * static_cast< double >( i ) );
        }
      }
    }

  typedef itk::Image< float, 4 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  mrfit::MRFit::Pointer fit = mrfit::MRFit::New();
  fit->SetInput( reader->GetOutput() );
  fit->SetTimes( echos );
  fit->FitT2( static_cast< mrfit::MRFit::T2FittingType >( algorithm ), maxT2 );

  typedef itk::Image< float, 3 > OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  WriterType::Pointer writer = WriterType::New();

  // output T2
  writer->SetInput( fit->GetMap( 0 ) );
  writer->SetFileName( outputFileName.c_str() );
  writer->Update();

  // output A
  if ( outputFileNameA != "" )
    {
    writer->SetInput( fit->GetMap( 1 ) );
    writer->SetFileName( outputFileNameA.c_str() );
    writer->Update();
    }

  // output C
  if ( outputFileNameC != "" )
    {
    writer->SetInput( fit->GetMap( 2 ) );
    writer->SetFileName( outputFileNameC.c_str() );
    writer->Update();
    }

  // output r^2
  if ( outputFileNameR2 != "" )
    {
    writer->SetInput( fit->GetMap( 3 ) );
    writer->SetFileName( outputFileNameR2.c_str() );
    writer->Update();
    }

  return 0;
}
