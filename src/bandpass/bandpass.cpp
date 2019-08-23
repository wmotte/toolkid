#include "tkdCmdParser.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template< class T >
void swap( T& a, T& b )
{
  T tmp = a;
  a = b;
  b = tmp;
}

/** (c) Press et al., Numerical Recipes 3rd ed., p.612 **/
template< class T >
void four1( T* data, const int& n, const int& isign )
{
  int nn, mmax, m, j, istep, i;
  T wtemp, wr, wpr, wpi, wi, theta, tempr, tempi, norm;

  if ( n < 2 || n & ( n - 1 ) )
    {
    throw( "n must be a power of 2" );
    }

  nn = n << 1;
  j = 1;

  for( i = 1; i < nn; i += 2 )
    {
    if ( j > i )
      {
      swap< T >( data[ j - 1 ], data[ i - 1 ] );
      swap< T >( data[ j ], data[ i ] );
      }
    m = n;
    while( m >= 2 && j > m )
      {
      j -= m;
      m >>= 1;
      }
    j += m;
    }

  mmax = 2;

  while( nn > mmax )
    {
    istep = mmax << 1;
    theta = static_cast< T >( isign ) * ( 6.28318530717959 / mmax );
    wtemp = sin( 0.5 * theta );
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin( theta );
    wr = 1.0;
    wi = 0.0;

    for( m = 1; m < mmax; m +=2 )
      {
      for( i = m; i <= nn; i += istep )
        {
        j = i + mmax;
        tempr = wr * data[ j - 1 ] - wi * data[ j ];
        tempi = wr * data[ j ] + wi * data[ j - 1 ];
        data[ j - 1 ] = data[ i - 1 ] - tempr;
        data[ j ] = data[ i ] - tempi;
        data[ i - 1 ] += tempr;
        data[ i ] += tempi;
        }
      wtemp = wr;
      wr = wtemp * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
      }
    mmax = istep;
    }

  if ( isign == -1 )
    {
    norm = 1. / n;
    for( i = 0; i < nn; ++i )
      {
      data[ i ] *= norm;
      }
    }
}

/** (c) Press et al., Numerical Recipes 3rd ed., p.784 **/
void fitline( const vnl_vector< float >& data, float& a, float& b )
{
  int numberOfPoints = data.size();

  float sx = 0;
  float sy = 0;
  float st2 = 0;

  for( int i = 0; i < numberOfPoints; ++i )
    {
    sx += static_cast< float >( i );
    sy += data( i );
    }

  float sxoss = sx / static_cast< float >( numberOfPoints );

  for( int i = 0; i < numberOfPoints; ++i )
    {
    float t = static_cast< float >( i ) - sxoss;
    st2 += t * t;
    a += t * data( i );
    }
  a /= st2;
  b = ( sy - sx * a ) / static_cast< float >( numberOfPoints );
}

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "bandpass", "Band-pass 4D image filter" );

  std::string inputFileName, outputFileName;
  float lowPass, highPass;
  float sampleSpacing = 0;
  bool retrend = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "4D input image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "4D filtered output image" )
    ->SetRequired( true );

  p.AddArgument( lowPass, "lowpass" )
    ->AddAlias( "lp" )
    ->SetInput( "float" )
    ->SetDescription( "Low pass (Hz)" )
    ->SetRequired( true );

  p.AddArgument( highPass, "highpass" )
    ->AddAlias( "hp" )
    ->SetInput( "float" )
    ->SetDescription( "High pass (Hz)" )
    ->SetRequired( true );

  p.AddArgument( sampleSpacing, "sample-spacing" )
    ->AddAlias( "ss" )
    ->SetInput( "float" )
    ->SetDescription( "Sample spacing (in seconds; default: 4th voxel dimension)" );

  p.AddArgument( retrend, "retrend" )
    ->AddAlias( "r" )
    ->SetInput( "boolean" )
    ->SetDescription( "Retrend data after filtering" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef itk::Image< float, 4 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  ImageType::Pointer image = reader->GetOutput();
  reader = 0;

  ImageType::RegionType region = image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::SpacingType spacing = image->GetSpacing();

  if ( sampleSpacing <= 0 )
    {
    sampleSpacing = spacing[ 3 ];
    }

  int numberOfVoxels = size[ 0 ] * size[ 1 ] * size[ 2 ];
  int numberOfTimePoints = size[ 3 ];

  int bufferLength = pow( 2, vcl_ceil( vcl_log( numberOfTimePoints ) / vcl_log( 2. ) ) );
  float sampleFactor = static_cast< float >( bufferLength ) * sampleSpacing;
  int lowerIndex1 = static_cast< int >( vcl_ceil( highPass * sampleFactor ) );
  int higherIndex1 = static_cast< int >( vcl_floor( lowPass * sampleFactor ) );
  int lowerIndex2 = bufferLength - higherIndex1 - 1;
  int higherIndex2 = bufferLength - lowerIndex1 - 1;

  vnl_matrix_ref< float > matrix(
      numberOfTimePoints,
      numberOfVoxels,
      image->GetPixelContainer()->GetBufferPointer() );

  vnl_vector< vcl_complex< float > > buffer( bufferLength );
  buffer.fill( itk::NumericTraits< vcl_complex< float > >::Zero );

  for( int i = 0; i < numberOfVoxels; ++i )
  {
    float a = 0;
    float b = 0;


    // FITLINE...
    fitline( matrix.get_column( i ), a, b );

    for( int j = 0; j < numberOfTimePoints; ++j )
      {
      float onLine = a * static_cast< float >( j ) + b;
      buffer( j ) = vcl_complex< float >( matrix( j, i ) - onLine, 0 );
      }

    // zero-filling to power of 2
    for( int j = numberOfTimePoints; j < bufferLength; ++j )
      {
      buffer( j ) = 0;
      }

    // forward transform
    four1< float >( reinterpret_cast< float* >( buffer.data_block() ), buffer.size(), 1 );

    for( int j = 0; j < bufferLength; ++j )
      {
      if ( j >= lowerIndex1 && j <= higherIndex1 )
        {
        continue;
        }

      if ( j >= lowerIndex2 && j <= higherIndex2 )
        {
        continue;
        }

      buffer( j ) = 0;
      }

    // backward transform
    four1< float >( reinterpret_cast< float* >( buffer.data_block() ), buffer.size(), -1 );

    for( int j = 0; j < numberOfTimePoints; ++j )
      {
      const vcl_complex< float >& p = buffer( j );
      matrix( j, i ) = p.real();

      if ( retrend )
        {
        matrix( j, i ) += a * static_cast< float >( j ) + b;
        }
      }
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( outputFileName.c_str() );
  writer->Update();

  return 0;
}
