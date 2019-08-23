#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "itkNumericTraits.h"
#include "tkdCmdParser.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkTimeProbe.h"

class RandomizeGraph
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 2 > ImageType;
  typedef vnl_matrix_ref< PixelType > DataMatrixType;
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  void Run( const std::string& filename, int numberOfIterations, const std::string& outputFileName, bool isUndirected )
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filename.c_str() );
    reader->Update();

    ImageType::Pointer image = reader->GetOutput();
    reader = 0;

    PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
    ImageType::RegionType region = image->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    int rows = size[ 0 ];
    int cols = size[ 1 ];

    DataMatrixType R( rows, cols, buffer );

    GeneratorType::Pointer generator = GeneratorType::New();
    generator->SetSeed();

    std::vector< int > indicesRow, indicesColumn;

    for( int j = 0; j < cols; ++j )
      {
      for( int i = 0; i < rows; ++i )
        {
        indicesRow.push_back( i );
        indicesColumn.push_back( j );

        if ( R( i, j ) < 0 )
          {
          R( i, j ) = 0;
          }
        }
      }

    std::vector< PixelType > diagonal = std::vector< PixelType >( rows );
    for( int i = 0; i < rows; ++i )
      {
      diagonal[ i ] = R( i, i );
      }

    int K = indicesRow.size();
    numberOfIterations *= K;

    std::cout << rows << " rows/columns" << std::endl;
    std::cout << "K=" << K << std::endl;
    std::cout << "Iterations: " << numberOfIterations << std::endl;

    int percent = numberOfIterations / 100;

    itk::TimeProbe probe;

    for( int i = 0; i < numberOfIterations; ++i )
      {
      probe.Start();

      if ( ( ( i + 1 ) % percent ) == 0 )
        {
        std::cout << ( i + 1 ) << "/" << numberOfIterations << ": " << probe.GetMeanTime() << "\r";
        std::cout.flush();
        }

      int a, b, c, d, e1, e2;

      while( true )
        {
        do
          {
          e1 = generator->GetIntegerVariate( K - 1 );

          do
            {
            e2 = generator->GetIntegerVariate( K - 1 );
            }
          while( e1 == e2 );

          a = indicesRow[ e1 ];
          b = indicesColumn[ e1 ];
          c = indicesRow[ e2 ];
          d = indicesColumn[ e2 ];

          }
        while ( a == c || a == d || b == c || b == d );

        if ( !( R( a, d ) || R( c, b ) ) )
          {
          R( a, d ) = R( a, b );
          R( a, b ) = 0;
          R( c, b ) = R( c, d );
          R( c, d ) = 0;

          indicesColumn[ e1 ] = d;
          indicesColumn[ e2 ] = b;
          break;
          }
        }
      probe.Stop();
      }

    std::cout << probe.GetMeanTime() << std::endl;

    if ( isUndirected )
      {
      // copy to lower triangle
      for( int i = 0; i < rows; ++i )
        {
        for( int j = i + 1; j < cols; ++j )
          {
          R( j, i ) = R( i, j );
          }
        }
      }

    for( int i = 0; i < rows; ++i )
      {
      R( i, i ) = diagonal[ i ];
      }

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputFileName.c_str() );
    writer->SetInput( image );
    writer->Update();
  }

  void RunDense( const std::string& filename, const std::string& outputFileName, bool isUndirected )
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filename.c_str() );
    reader->Update();

    ImageType::Pointer image = reader->GetOutput();
    reader = 0;

    PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
    ImageType::RegionType region = image->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    int rows = size[ 0 ];
    int cols = size[ 1 ];

    DataMatrixType R( rows, cols, buffer );

    GeneratorType::Pointer generator = GeneratorType::New();
    generator->SetSeed();

    int numberOfValues = rows * cols;
    std::vector< PixelType > values = std::vector< PixelType >( numberOfValues );
    for( int i = 0; i < numberOfValues; ++i )
      {
      values[ i ] = buffer[ i ];
      }

    std::vector< PixelType > diagonal = std::vector< PixelType >( rows );
    for( int i = 0; i < rows; ++i )
      {
      diagonal[ i ] = R( i, i );
      }

    for( int i = 0; i < numberOfValues; ++i )
      {
      int index = generator->GetIntegerVariate( numberOfValues - i - 1 );
      buffer[ i ] = values[ index ];
      values[ index ] = values[ numberOfValues - i - 1 ];
      }

    if ( isUndirected )
      {
      // copy to lower triangle
      for( int i = 0; i < rows; ++i )
        {
        for( int j = i + 1; j < cols; ++j )
          {
          R( j, i ) = R( i, j );
          }
        }
      }

    for( int i = 0; i < rows; ++i )
      {
      R( i, i ) = diagonal[ i ];
      }

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputFileName.c_str() );
    writer->SetInput( image );
    writer->Update();
  }
};

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "randomizegraph", "Randomize weighted graphs" );

  std::string inputFileName;
  std::string outputFileName;
  int numberOfIterations = 10;
  bool isDenseGraph = false;
  bool isUndirected = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input image: 2D adjacency matrix" )
    ->SetRequired( true );

  p.AddArgument( numberOfIterations, "iterations" )
    ->AddAlias( "n" )
    ->SetDescription( "Number of iterations for sparse graphs (default: 10)" );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output image: 2D adjacency matrix" )
    ->SetRequired( true );

  p.AddArgument( isDenseGraph, "dense" )
    ->AddAlias( "d" )
    ->SetDescription( "Input graph is fully connected with weights; randomize by weights shuffling" );

  p.AddArgument( isUndirected, "undirected" )
    ->AddAlias( "u" )
    ->SetDescription( "Keep graph undirected by symmetrizing the final adjacency matrix" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  RandomizeGraph randomize;

  if ( isDenseGraph )
    {
    randomize.RunDense( inputFileName, outputFileName, isUndirected );
    }
  else
    {
    randomize.Run( inputFileName, numberOfIterations, outputFileName, isUndirected );
    }

  return 0;
}
