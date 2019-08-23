#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_vector.h"
#include "itkImageFileReader.h"
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "meants", "Mean signal from 4D time series image" );

  std::string inputImageName;
  std::string maskImageName;
  std::string outputFileName;
  std::string labelsFileName;
  std::string xLabelsFileName;
  std::vector< int > selectedLabelsList;
  double spacingX = 1;
  double startX = 0;
  bool noHeader = false;
  bool noX = false;
  bool noSD = false;
  bool log = false;
  bool sub1 = false;

  p.AddArgument( inputImageName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input 4D time series image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output matrix (CSV)" )
    ->SetRequired( true );

  p.AddArgument( maskImageName, "mask" )
    ->AddAlias( "m" )
    ->SetDescription( "Input mask/label image" )
    ->SetRequired( true );

  p.AddArgument( labelsFileName, "labels" )
    ->AddAlias( "l" )
    ->SetDescription( "Text file containing label names" );

  p.AddArgument( spacingX, "x-spacing" )
    ->AddAlias( "sp" )
    ->SetDescription( "x-axis spacing (default: 1)" );

  p.AddArgument( startX, "x-start" )
    ->AddAlias( "st" )
    ->SetDescription( "x-axis start value (default: 0)" );

  p.AddArgument( selectedLabelsList, "select" )
    ->SetDescription( "Select a subset of labels" );

  p.AddArgument( xLabelsFileName, "x-labels" )
    ->AddAlias( "xl" )
    ->SetDescription( "x-axis labels file name" );

  p.AddArgument( noHeader, "no-header" )
    ->AddAlias( "nh" )
    ->SetDescription( "Do not output header" );

  p.AddArgument( noX, "no-x-axis" )
    ->AddAlias( "nx" )
    ->SetDescription( "Do not output x-axis values" );

  p.AddArgument( noSD, "no-sd" )
    ->AddAlias( "nsd" )
    ->SetDescription( "Do not output standard deviations" );

  p.AddArgument( log, "logarithm" )
    ->AddAlias( "log" )
    ->SetDescription( "Calculate mean from log-transformed values" );

  p.AddArgument( sub1, "subtract-one" )
    ->AddAlias( "s1" )
    ->SetDescription( "Subtract 1 from input" );

  if ( !p.Parse( argc, argv ) )
  {
    p.PrintUsage( std::cout );
    return -1;
  }

  std::set< int > selectedLabels;
  for( std::vector< int >::iterator i = selectedLabelsList.begin(); i != selectedLabelsList.end(); ++i )
  {
    selectedLabels.insert( *i );
  }

  typedef itk::Image< float, 4 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::Image< int, 3 > LabelType;
  typedef itk::ImageFileReader< LabelType > LabelReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageName.c_str() );
  reader->Update();

  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( maskImageName.c_str() );
  labelReader->Update();

  std::map< int, std::string > labelNames;

  if ( labelsFileName != "" )
  {
    std::ifstream in( labelsFileName.c_str() );
    while( !in.eof() )
    {
      int key;
      std::string name;

      std::string buffer;
      std::getline( in, buffer, '\t' );
      if ( !in.good() )
      {
        break;
      }

      std::stringstream ss;
      ss << buffer;
      ss >> key;

      std::getline( in, name );
      labelNames[ key ] = name;
    }
  }


  ImageType::Pointer image = reader->GetOutput();
  LabelType::Pointer label = labelReader->GetOutput();

  reader = 0;
  labelReader = 0;

  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  int numberOfTimepoints = size[ 3 ];
  int numberOfVoxels = size[ 0 ] * size[ 1 ] * size[ 2 ];

  std::vector< double > xLabels;
  if ( xLabelsFileName != "" )
  {
    std::ifstream in( xLabelsFileName.c_str() );
    while( !in.eof() )
    {
      double x;
      in >> x;
      if ( in.good() )
      {
        xLabels.push_back( x );
      }
    }
  }
  else
  {
    for( int i = 0; i < numberOfTimepoints; ++i )
    {
      xLabels.push_back( startX + static_cast< double >( i ) * spacingX );
    }
  }

  vnl_matrix_ref< float > matrix( numberOfTimepoints, numberOfVoxels, image->GetPixelContainer()->GetBufferPointer() );
  int* mask = label->GetPixelContainer()->GetBufferPointer();

  if ( log )
  {
    for( int i = 0; i < numberOfVoxels; ++i )
    {
      for( int j = 0; j < numberOfTimepoints; ++j )
      {
        matrix( j, i ) = vcl_log( matrix( j, i ) );
      }
    }
  }

  std::map< int, int > labels;
  int numberOfLabels = 0;
  for( int i = 0; i < numberOfVoxels; ++i )
  {
    if ( mask[ i ] == 0 || labels.find( mask[ i ] ) != labels.end() )
    {
      continue;
    }

    if ( selectedLabels.size() > 0 && selectedLabels.find( mask[ i ] ) != selectedLabels.end() )
    {
      mask[ i ] = 0;
      continue;
    }

    labels[ mask[ i ] ] = numberOfLabels++;
  }

  vnl_matrix< float > means( numberOfTimepoints, numberOfLabels );
  means.fill( 0 );

  vnl_vector< int > count( numberOfLabels );
  count.fill( 0 );

  for( int i = 0; i < numberOfVoxels; ++i )
  {
    if ( mask[ i ] == 0 )
    {
      continue;
    }

    int column = labels[ mask[ i ] ];
    vnl_vector< float > signal = matrix.get_column( i );
    means.set_column( column, means.get_column( column ) + signal );
    count[ column ]++;
  }

  for( int i = 0; i < numberOfLabels; ++i )
  {
    means.set_column( i, means.get_column( i ) / static_cast< float >( count[ i ] ) );
  }

  vnl_matrix< float > sd( numberOfTimepoints, numberOfLabels );
  sd.fill( 0 );

  for( int i = 0; i < numberOfVoxels; ++i )
  {
    if ( mask[ i ] == 0 )
    {
      continue;
    }

    int column = labels[ mask[ i ] ];
    vnl_vector< float > signal = matrix.get_column( i ) - means.get_column( column );
    for( int j = 0; j < numberOfTimepoints; ++j )
    {
      sd( j, column ) += ( signal( j ) * signal( j ) );
    }
  }

  for( int i = 0; i < numberOfTimepoints; ++i )
  {
    for( int j = 0; j < numberOfLabels; ++j )
    {
      sd( i, j ) = vcl_sqrt( sd( i, j ) / static_cast< float >( count[ j ] - 1 ) );
    }
  }

  std::ofstream out( outputFileName.c_str() );

  if ( labelNames.size() > 0 && !noHeader )
  {
    for( std::map< int, int >::iterator i = labels.begin(); i != labels.end(); ++i )
    {
      out << ( i == labels.begin() && noX ? "" : "," ) << "\"mean(" << labelNames[ i->first ] << ")\"";
      if ( !noSD )
      {
        out << ",\"sd(" << labelNames[ i->first ] << ")\"";
      }
    }
  }

  for( int i = 0; i < numberOfTimepoints; ++i )
  {
    if ( i > 0 || !noHeader )
    {
      out << std::endl;
    }

    if ( !noX )
    {
      out << xLabels[ i ];
    }

    for( int j = 0; j < numberOfLabels; ++j )
    {
      out << ( noX && j == 0 ? "" : "," ) << ( ( log ? vcl_exp( means( i, j ) ) - 1. : means( i, j ) ) - ( sub1 ? 1 : 0 ) );
      if ( !noSD )
      {
        out << "," << ( log ? vcl_exp( sd( i, j ) ) - 1. : sd( i, j ) );
      }
    }
  }

  out << std::endl;

  return 0;
}

