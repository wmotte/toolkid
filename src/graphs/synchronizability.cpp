#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "itkNumericTraits.h"
#include "tkdCmdParser.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

class GraphSynchronizability
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 2 > ImageType;
  typedef vnl_matrix_ref< PixelType > DataMatrixType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef vnl_symmetric_eigensystem< PixelType > EigenSystem;

  void CalculateLaplacian( DataMatrixType& G )
  {
    int rows = G.rows();
    int cols = G.cols();

    // calculate node degrees
    std::vector< PixelType > degrees( rows );
    for( int i = 0; i < rows; ++i )
      {
      PixelType d = 0;
      for( int j = 0; j < cols; ++j )
        {
        d += G( i, j );
        }
      degrees[ i ] = d;
      }

    // D: diagonal degrees matrix
    // A: adjacency matrix
    // L: Laplacian matrix, such that:
    // L = D - A
    for( int i = 0; i < rows; ++i )
      {
      for( int j = 0; j < cols; ++j )
        {
        G( i, j ) *= -1.;
        }

      G( i, i ) += degrees[ i ];
      }
  }

  PixelType CalculateSynchronizability( DataMatrixType& G )
  {
    EigenSystem es( G );

    // ratio of largest eigenvalue and second smallest eigenvalue (since smallest ~ 0)
    return es.get_eigenvalue( G.rows() - 1 ) / es.get_eigenvalue( 1 );
  }

  PixelType Run( const std::string& filename )
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

    DataMatrixType G( rows, cols, buffer );
    this->CalculateLaplacian( G );

    return this->CalculateSynchronizability( G );
  }

};

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "graphsynchronizability", "Compute synchronizability value R for a symmetric adjacency matrix" );

  std::string inputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input image: 2D adjacency matrix" )
    ->SetRequired( true );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  GraphSynchronizability gs;
  std::cout << "Synchronizability," << gs.Run( inputFileName ) << std::endl;

  return 0;
}
