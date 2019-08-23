#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "itkImageFileReader.h"
#include <iostream>
#include <fstream>
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "matrix2csv", "Export 2D matrix image to spreadsheet (CSV)" );

  std::string inputFileName, outputFileName, labelsFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input image: 2D matrix" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output csv spreadsheet" )
    ->SetRequired( true );

  p.AddArgument( labelsFileName, "labels" )
    ->AddAlias( "l" )
    ->SetInput( "filename" )
    ->SetDescription( "Input labels file" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef float PixelType;
  typedef itk::Image< PixelType, 2 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef vnl_matrix_ref< PixelType > MatrixType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  ImageType::Pointer image = reader->GetOutput();
  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();

  PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
  MatrixType matrix( size[ 0 ], size[ 1 ], buffer );

  std::vector< std::string > labels;
  if ( labelsFileName != "" )
    {
    std::ifstream in( labelsFileName.c_str() );
    while( !in.eof() )
      {
      std::string line;
      std::getline( in, line );

      if ( line.length() > 0 && in.good() && !in.fail() && !in.bad() )
        {
        labels.push_back( line );
        }
      }
    }

  std::ofstream out( outputFileName.c_str() );
  for( int i = 0; i < labels.size(); ++i )
    {
    out << ",\"" << labels[ i ] << "\"";
    }

  if ( labels.size() > 0 )
    {
    out << std::endl;
    }

  for( int i = 0; i < size[ 0 ]; ++i )
    {
    if ( i < labels.size() )
      {
      out << labels[ i ] << ( size[ 1 ] > 0 ? "," : "" );
      }

    for( int j = 0; j < size[ 1 ]; ++j )
      {
      out << ( j > 0 ? "," : "" );
      if ( matrix( i, j ) != matrix( i, j ) )
        {
        continue;
        }

      out << matrix( i, j );
      }
    out << std::endl;
    }

  return 0;
}

