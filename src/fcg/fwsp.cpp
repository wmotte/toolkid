#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "itkNumericTraits.h"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/floyd_warshall_shortest.hpp"
#include "boost/graph/adjacency_matrix.hpp"
#include "itkTimeProbe.h"
#include "tkdCmdParser.h"

template< class DirectedProperty=boost::undirectedS >
class ShortestPath
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 2 > ImageType;
  typedef vnl_matrix_ref< PixelType > DataMatrixType;

  void Run( const std::string& filename, PixelType threshold, const std::string& outputFileName )
  {
    typedef itk::ImageFileReader< ImageType > ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filename.c_str() );
    reader->Update();

    typename ImageType::Pointer image = reader->GetOutput();
    reader = 0;

    PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
    typename ImageType::RegionType region = image->GetLargestPossibleRegion();
    typename ImageType::SizeType size = region.GetSize();
    int rows = size[ 0 ];
    int cols = size[ 1 ];

    DataMatrixType data( rows, cols, buffer );

    typename ImageType::Pointer output = ImageType::New();
    output->CopyInformation( image );
    output->SetRegions( image->GetLargestPossibleRegion() );
    output->Allocate();
    output->FillBuffer( 0 );

    DataMatrixType distance( rows, cols, output->GetPixelContainer()->GetBufferPointer() );

    typedef boost::property< boost::edge_weight_t, PixelType > EdgeProperty;
    typedef boost::no_property VertexProperty;
    typedef boost::adjacency_matrix< DirectedProperty, VertexProperty, EdgeProperty > GraphType;

    GraphType graph( rows );

    typename boost::property_map< GraphType, boost::edge_weight_t >::type weightmap = boost::get( boost::edge_weight, graph );

    const bool isUndirected = typeid( DirectedProperty ) == typeid( boost::undirectedS );
    for( int i = 0; i < rows; ++i )
      {
      int start = isUndirected ? i + 1 : 0;
      for( int j = start; j < cols; ++j )
        {
        if ( i == j )
          {
          continue;
          }

        // absolute value of correlation coefficient in [0, 1]
        PixelType weight = vnl_math_abs( data( i, j ) );

        if ( weight <= threshold && weight > 0 )
          {
          typename GraphType::edge_descriptor e;
          bool inserted;
          boost::tie( e, inserted ) = boost::add_edge( i, j, graph );

          // insert inverse correlation coefficient as weight
          weightmap[ e ] = 1.0 / weight;
          }
        }
      }

    itk::TimeProbe probe;
    probe.Start();

    boost::floyd_warshall_all_pairs_shortest_paths(
        graph,
        distance );

    probe.Stop();

    PixelType sum = itk::NumericTraits< PixelType >::Zero;
    PixelType infinity = itk::NumericTraits< PixelType >::max();
    int count = 0;
    for( int i = 0; i < rows; ++i )
      {
      for( int j = 0; j < cols; ++j )
        {
        if ( i == j || distance( i, j ) == infinity || distance( i, j ) != distance( i, j ) )
          {
          continue;
          }

        ++count;
        sum += ( 1. / distance( i, j ) );
        }
      }

    PixelType mean = static_cast< PixelType >( count  ) / sum;

    typedef itk::ImageFileWriter< ImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputFileName.c_str() );
    writer->SetInput( output );
    writer->Update();

    std::cout << mean << std::endl;
  }
};

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "fwsp", "Floyd-Warshall weighted shortest paths on full graph" );

  std::string inputFileName;
  std::string outputFileName;
  float threshold = 1.0;
  bool isUndirected = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input image: 2D adjacency matrix with correlation coefficients" )
    ->SetRequired( true );

  p.AddArgument( threshold, "threshold" )
    ->AddAlias( "t" )
    ->SetInput( "float" )
    ->SetDescription( "Threshold; only include paths with weight in (0, threshold] (default: 1.0)" );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output image: 2D shortest paths matrix" )
    ->SetRequired( true );

  p.AddArgument( isUndirected, "undirected" )
    ->AddAlias( "u" )
    ->SetDescription( "Input is a symmetric matrix representing an undirected graph" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  if ( isUndirected )
    {
    ShortestPath< boost::undirectedS > sp;
    sp.Run( inputFileName, threshold, outputFileName );
    }
  else
    {
    ShortestPath< boost::directedS > sp;
    sp.Run( inputFileName, threshold, outputFileName );
    }

  return 0;
}
