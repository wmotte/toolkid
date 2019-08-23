#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"
#include "itkPermuteAxesImageFilter.h"
#include "tkdCmdParser.h"

#define logMacro(x) \
  { \
  if ( m_Verbose ) \
    { \
      std::cout << x << std::endl; \
    } \
  }

class ExportBFloat
{
private:

bool m_Verbose;
std::string m_InputFileName;
std::string m_OutputFileName;

void ByteSwap( char* block, int size )
{
  register int i = 0;
  register int j = size - 1;

  for( ; i < j ; ++i, --j )
  {
    std::swap( block[ i ], block[ j ] );
  }
}

public:

int Run( int argc, char ** argv )
{
  std::map< std::string, std::vector< std::string > > args;

  tkd::CmdParser p( "nii2bfloat", "Convert nifti image to bfloat(s)" );
  std::string inputFileName;
  std::string outputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename(s)" )
    ->SetDescription( "Input image(s)" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output image" )
    ->SetRequired( true );

  p.AddArgument( m_Verbose, "verbose" )
    ->AddAlias( "v" )
    ->SetDescription( "Verbose" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  if ( outputFileName.length() > 7 && outputFileName.substr( outputFileName.length() - 6, 6 ) == "bfloat" )
    {
    outputFileName = outputFileName.substr( 0, outputFileName.length() - 7 );
    }

  typedef itk::Image< float, 4 > ImageType;
  typedef itk::Image< float, 4 > SliceType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ExtractImageFilter< ImageType, SliceType > ExtractType;
  typedef itk::PermuteAxesImageFilter< ImageType > FilterType;

  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();

  reader->SetFileName( inputFileName.c_str() );
  filter->SetInput( reader->GetOutput() );

  reader->Update();
  int volumes = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[ 3 ];

  FilterType::PermuteOrderArrayType order;

  if ( volumes > 1 )
    {
    order[ 0 ] = 0;
    order[ 1 ] = 1;
    order[ 2 ] = 3;
    order[ 3 ] = 2;
    }
  else
    {
    order[ 0 ] = 0;
    order[ 1 ] = 1;
    order[ 2 ] = 2;
    order[ 3 ] = 3;
    }

  filter->SetOrder( order );

  filter->Update();

  ImageType::Pointer image = filter->GetOutput();
  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();

  ImageType::RegionType selectRegion;
  ImageType::SizeType selectSize;
  ImageType::IndexType selectIndex;

  selectSize[ 0 ] = size[ 0 ];
  selectSize[ 1 ] = size[ 1 ];
  selectSize[ 2 ] = size[ 2 ];
  selectSize[ 3 ] = 1;

  selectIndex[ 0 ] = 0;
  selectIndex[ 1 ] = 0;
  selectIndex[ 2 ] = 0;

  int lastLength = 0;
  for( int n = size[ 3 ] - 1; n >= 0; --n )
    {
    ExtractType::Pointer extract = ExtractType::New();
    extract->SetInput( filter->GetOutput() );

    selectIndex[ 3 ] = n;

    selectRegion.SetSize( selectSize );
    selectRegion.SetIndex( selectIndex );

    extract->SetExtractionRegion( selectRegion );
    extract->Update();

    float* data = extract->GetOutput()->GetPixelContainer()->GetBufferPointer();
    unsigned int size = extract->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels() * 4;
    int bufferSize = size / 4;

    logMacro( "Byte swap block of size " << size << " bytes" );

    for( int i = 0; i < bufferSize; ++i )
      {
      ByteSwap( reinterpret_cast< char* >( & ( data[ i ] ) ), 4 );
      }

    std::stringstream numberStream;
    std::stringstream preStream;

    numberStream << n;

    int length = numberStream.str().length();
    if ( length < lastLength )
      {
      int difference = lastLength - length;
      for( int j = 0; j < difference; ++j )
        {
        preStream << "0";
        }
      }

    std::stringstream fileNumber;
    fileNumber << preStream.str() << numberStream.str();
    lastLength = fileNumber.str().length();

    std::stringstream outputFileNameHeader, outputFileNameBlock;

    if ( volumes > 1 )
      {
      outputFileNameHeader << outputFileName << "_" << fileNumber.str() << ".hdr";
      outputFileNameBlock << outputFileName << "_" << fileNumber.str() << ".bfloat";
      }
    else
      {
      outputFileNameHeader << outputFileName << ".hdr";
      outputFileNameBlock << outputFileName << ".bfloat";
      }

    std::ofstream out( outputFileNameBlock.str().c_str(), std::ios_base::binary );
    out.write( reinterpret_cast< char* >( data ), size );

    SliceType::RegionType sliceRegion = extract->GetOutput()->GetLargestPossibleRegion();
    SliceType::SizeType sliceSize = sliceRegion.GetSize();

    std::ofstream outHeader( outputFileNameHeader.str().c_str() );
    outHeader << sliceSize[ 1 ] << " " << sliceSize[ 0 ] << " " << sliceSize[ 2 ];
    outHeader << std::endl;
  }

  return 0;
}


ExportBFloat()
{
  m_Verbose = false;
}

};

int main( int argc, char ** argv )
{
  ExportBFloat e;

  return e.Run( argc, argv );
}

