#include "itkImage.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix_ref.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "tkdCmdParser.h"
#include "nrFourier.h"
#include "nrCorrelation.h"
#include <sstream>

template< class T >
void Hilbert( const vnl_vector< T >& input, const int& n, vnl_vector< T >& output )
{
  const int n2 = 2 * n;
  output = input;
  nr::fourier::FFT< T >::four1( output.data_block(), n, 1 );

  // k=1 ... n/2
  for( int k = 2; k < n; ++k )
    {
    output( k ) *= 2.;
    }
    
  // k = n/2+1 ... n  
  for( int k = n + 2; k < n2; ++k )
    {
    output( k ) = 0;
    }
    
  nr::fourier::FFT< T >::four1( output.data_block(), n, -1 );
  
  for( int i = 1; i < n2; i += 2 )
    {
    output( i ) *= -1.;
    }
}

template< class T >
void Phase( const vnl_vector< T >& input, vnl_vector< T >& output )
{
  const int n = input.size() / 2;
  
  output.set_size( n );
  for( int i = 0; i < n; ++i )
    {
    output( i ) = atan2( input[ i * 2 + 1 ], input[ i * 2 ] );
    }
}

template< class T >
T PhaseLagIndex( const vnl_vector< T >& phase1, const vnl_vector< T >& phase2 )
{
  const int n = phase1.size();
  
  int count1 = 0;
  int count2 = 0;
  for( int i = 0; i < n; ++i )
    {
    const T difference = phase1( i ) - phase2( i );
    
    if ( difference > 0 )
      {
      ++count1;
      ++count2;
      }
    else if ( difference < 0 )
      {
      ++count2;
      }
    }
    
  if ( count2 == 0 )
    {
    return 0;
    }
    
  return 2. * vnl_math_abs( 0.5 - static_cast< T >( count1 ) / static_cast< T >( count2 ) );
}

template< class T >
void RealToComplexVector( const vnl_vector< T >& input, vnl_vector< T >& output, int n )
{
  output.set_size( n * 2 );
  output.fill( 0 );
  
  for( unsigned int i = 0; i < input.size(); ++i )
    {
    output[ i * 2 ] = input[ i ];
    }
}

template< class T >
T PLI( const vnl_vector< T >& a, const vnl_vector< T >& b, int n )
{
  vnl_vector< T > vectorA, vectorB;
  vnl_vector< T > hilbertA, hilbertB;
  vnl_vector< T > phaseA, phaseB;

  RealToComplexVector< T >( a, vectorA, n );
  RealToComplexVector< T >( b, vectorB, n );
  
  Hilbert< T >( vectorA, n, hilbertA );
  Hilbert< T >( vectorB, n, hilbertB );

  Phase< T >( hilbertA, phaseA );
  Phase< T >( hilbertB, phaseB );
  
  return PhaseLagIndex< T >( phaseA, phaseB );  
}


namespace fcg
{

class PLI
{
public:
  typedef float PixelType;
  typedef float OutputPixelType;
  typedef unsigned char MaskPixelType;

  typedef itk::Image< PixelType, 4 > ImageType;
  typedef itk::Image< MaskPixelType, 3 > MaskImageType;
  typedef itk::ImageMaskSpatialObject< 3 > MaskType;

  typedef vnl_vector< PixelType > VectorType;
  typedef vnl_matrix_ref< PixelType > DataType;
  typedef vnl_matrix< PixelType > MatrixType;
  typedef itk::Image< PixelType, 3 > MapType;
  typedef vnl_matrix_ref< MaskPixelType > DataMaskType;
  typedef itk::Image< OutputPixelType, 2 > OutputImageType;
  typedef vnl_matrix_ref< OutputPixelType > OutputDataType;

void Run( const std::string& inputFileName, const std::string& maskFileName, const std::string& outputFileName, int skipBegin = 50, int skipEnd = 50 )
{
  std::cout << "Read " << inputFileName << std::endl;
  ImageType::Pointer image = LoadImage( inputFileName );
  std::cout << "Read " << maskFileName << std::endl;
  MaskType::Pointer mask = LoadMask( maskFileName );

  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();

  int voxels = size[ 0 ] * size[ 1 ] * size[ 2 ];
  int timepoints = size[ 3 ];

  PixelType* data = image->GetPixelContainer()->GetBufferPointer();
  MaskPixelType* dataMask = const_cast< MaskType::ImageType* >( mask->GetImage() )->GetPixelContainer()->GetBufferPointer();

  DataType source( timepoints, voxels, data );
  MatrixType matrix = source.extract( source.rows() - skipBegin - skipEnd, source.cols(), skipBegin, 0 );

  std::vector< int > indices;
  std::vector< int > inverseIndices;
  int rowCount = 0;
  for( int i = 0; i < voxels; ++i )
    {
    if ( dataMask[ i ] )
      {
      indices.push_back( i );
      inverseIndices.push_back( rowCount++ );
      }
    else
      {
      inverseIndices.push_back( -1 );
      }
    }

  int validVoxels = indices.size();
  timepoints = matrix.rows();
  std::cout << validVoxels << " voxels with " << matrix.rows() << " timepoints" << std::endl;

  OutputImageType::Pointer outputImage = OutputImageType::New();
  OutputImageType::RegionType outputRegion;
  OutputImageType::SizeType outputSize;
  outputSize[ 0 ] = validVoxels;
  outputSize[ 1 ] = validVoxels;
  outputRegion.SetSize( outputSize );
  outputImage->SetRegions( outputRegion );
  outputImage->Allocate();
  outputImage->FillBuffer( itk::NumericTraits< OutputPixelType >::Zero );
  OutputPixelType* outputData = outputImage->GetPixelContainer()->GetBufferPointer();
  OutputDataType output( validVoxels, validVoxels, outputData );

  int n = pow( 2, vcl_ceil( vcl_log( matrix.rows() ) / vcl_log( 2. ) ) );
  int n2 = n * 2;

  
  int numberOfCalculations = ( validVoxels - 1 ) * validVoxels * 0.5 + validVoxels;
  int counter = 0;
  int lastPercentage = 0;
  
  MatrixType phase( n, validVoxels );

  vnl_vector< PixelType > inputVector, hilbertVector, phaseVector;
  
  for( int i = 0; i < validVoxels; ++i )
    {
    RealToComplexVector< PixelType >( matrix.get_column( indices[ i ] ), inputVector, n );
    Hilbert< PixelType >( inputVector, n, hilbertVector );
    Phase< PixelType >( hilbertVector, phaseVector );
    phase.set_column( i, phaseVector );
    
/*
    VectorType vector( n2 );
    vector.fill( 0 );
    const VectorType& a = matrix.get_column( indices[ i ] );
    
    // fill buffer with n x real,imag ... real,imag pairs
    PixelType* buffer = vector.data_block();
    for( int k = 0; k < timepoints; ++k )
      {
      *( buffer++ ) = a( k );
      *( buffer++ ) = 0;
      }
      
    // forward FFT
    nr::fourier::FFT< float >::four1( vector.data_block(), n, 1 );

    // zero negative frequencies

    // k=1 ... n/2
    for( int k = 2; k < n; ++k )
      {
      vector( k ) *= 2.;
      }
    
    // k = n/2+1 ... n  
    for( int k = n + 2; k < n2; ++k )
      {
      vector( k ) = 0;
      }
    
    // back-transform (yields conjugate)
    nr::fourier::FFT< float >::four1( vector.data_block(), n, -1 );
    
    // calculate instantaneous phases
    buffer = vector.data_block();
    for( int k = 0; k < n; ++k )
      {
      phase( k, i ) = vcl_atan2( - buffer[ 1 ], buffer[ 0 ] );
      buffer += 2;
      }
*/
    int percentage = ( static_cast< double >( ++counter ) / static_cast< double >( numberOfCalculations ) ) * 100.;
    if ( percentage > lastPercentage )
      {
      std::cout << "\r" << percentage << "%";
      std::cout.flush();
      lastPercentage = percentage;
      }
    }
  
  for( int i = 0; i < validVoxels; ++i )
    {
    const VectorType& phase1 = phase.get_column( i );
    
    for( int j = i + 1; j < validVoxels; ++j )
      {
      const VectorType& phase2 = phase.get_column( j );
      output( i, j ) = PhaseLagIndex< PixelType >( phase1, phase2 );
      output( j, i ) = output( i, j );
/*      
      VectorType difference = phase1 - phase2;
      int count1 = 0;
      int count2 = 0;
      
      for( int k = 0; k < n; ++k )
        {
        if ( difference( k ) > 0 )
          {
          ++count1;
          ++count2;
          }
        else if ( difference( k ) < 0 )
          {
          ++count2;
          }
        }
        
      output( i, j ) = count2 > 0 ? vnl_math_abs( 0.5 - static_cast< PixelType >( count1 ) / static_cast< PixelType >( count2 ) ) * 2. : 0;
      output( j, i ) = output( i, j );
*/      
      int percentage = ( static_cast< double >( ++counter ) / static_cast< double >( numberOfCalculations ) ) * 100.;
      if ( percentage > lastPercentage )
        {
        std::cout << "\r" << percentage << "%";
        std::cout.flush();
        lastPercentage = percentage;
        }
      }
    }
    
  std::cout << std::endl;

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( outputImage );

  std::cout << "Write " << outputFileName << std::endl;
  writer->Update();
}

ImageType::Pointer LoadImage( const std::string& filename )
{
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();

  return reader->GetOutput();
}

MaskType::Pointer LoadMask( const std::string& filename )
{
  typedef itk::ImageFileReader< MaskImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();

  MaskType::Pointer mask = MaskType::New();
  mask->SetImage( reader->GetOutput() );
  return mask;
}

};

} // end namespace fcg


int main( int argc, char ** argv )
{
  tkd::CmdParser p( "rmatrix", "Calculate matrix of voxel-wise correlation coefficients" );

  std::string inputFileName;
  std::string maskFileName;
  std::string outputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input 4D time series image" )
    ->SetRequired( true );

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetInput( "filename" )
    ->SetDescription( "Input 3D mask image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output 2D image: matrix of correlation coefficients" )
    ->SetRequired( true );

  std::vector< int > skip;
  skip.push_back( 50 );
  skip.push_back( 50 );

  p.AddArgument( skip, "skip" )
    ->AddAlias( "s" )
    ->SetMinMax( 2, 2 )
    ->SetDescription( "Skip samples at begin/end (default: 50 50)" );
    
  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  fcg::PLI pli;
  pli.Run( inputFileName, maskFileName, outputFileName, skip[ 0 ], skip[ 1 ] );

  return 0;
}
