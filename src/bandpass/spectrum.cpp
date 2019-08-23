#include "tkdCmdParser.h"
#include "vnl/algo/vnl_fft_1d.h"
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
void four1( vnl_vector< T >& data, const int& isign )
{
  int nn, mmax, m, j, istep, i;
  int n = data.size() / 2;
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

/** (c) Press et al., Numerical Recipes 3rd ed., p.654 **/
template< class T >
void realft( vnl_vector< T >& data, const int& isign )
{
  int i, i1, i2, i3, i4;
  int n = data.size();
  T c1 = 0.5;
  T c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp;
  T theta = 3.141592653589793238 / T( n >> 1 );
  if ( isign == 1 )
    {
    c2 = -0.5;
    four1< T >( data, 1 );
    }
  else
    {
    c2 = 0.5;
    theta = -theta;
    }

  wtemp = sin( 0.5 * theta );
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin( theta );
  wr = 1.0 + wpr;
  wi = wpi;

  for( i = 1; i < ( n >> 2 ); ++i )
    {
    i1 = i + i;
    i2 = 1 + i1;
    i3 = n - i1;
    i4 = 1 + i3;
    h1r = c1 * ( data[ i1 ] + data[ i3 ] );
    h1i = c1 * ( data[ i2 ] - data[ i4 ] );
    h2r = -c2 * ( data[ i2 ] + data[ i4 ] );
    h2i = c2 * ( data[ i1 ] - data[ i3 ] );
    data[ i1 ] = h1r + wr * h2r - wi * h2i;
    data[ i2 ] = h1i + wr * h2i + wi * h2r;
    data[ i3 ] = h1r - wr * h2r + wi * h2i;
    data[ i4 ] = -h1i + wr * h2i + wi * h2r;
    wtemp = wr;
    wr = wtemp * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
    }

  if ( isign == 1 )
    {
    h1r = data[ 0 ];
    data[ 0 ] = h1r + data[ 1 ];
    data[ 1 ] = h1r - data[ 1 ];
    }
  else
    {
    h1r = data[ 0 ];
    data[ 0 ] = c1 * ( h1r + data[ 1 ] );
    data[ 1 ] = c1 * ( h1r - data[ 1 ] );
    four1< T >( data, -1 );
    }
}

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "spectrum", "Spectrum from 4D image" );

  std::string inputFileName, outputFileName;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "4D input image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "4D magnitude output image" );

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

  int numberOfVoxels = size[ 0 ] * size[ 1 ] * size[ 2 ];
  int numberOfTimePoints = size[ 3 ];
  int bufferLength = pow( 2., vcl_ceil( vcl_log( numberOfTimePoints ) / vcl_log( 2. ) ) );

  vnl_matrix_ref< float > matrix(
      numberOfTimePoints,
      numberOfVoxels,
      image->GetPixelContainer()->GetBufferPointer() );

  vnl_vector< float > data( bufferLength );
  vnl_vector< float > output( bufferLength / 2 );

  for( int i = 0; i < numberOfVoxels; ++i )
    {
    for( int j = 0; j < numberOfTimePoints; ++j )
      {
      data( j ) = matrix( j, i );
      }

    for( int j = numberOfTimePoints; j < bufferLength; ++j )
      {
      data( j ) = 0;
      }

    realft< float >( data, 1 );
    float fac = 2. / static_cast< float >( bufferLength * bufferLength );

    output( 0 ) = 0.5 * fac * data[ 0 ] * data[ 0 ];
    for( int j = 1; j < bufferLength / 2; ++j )
      {
      output( j ) = fac * ( data[ 2 * j ] * data[ 2 * j ] + data[ 2 * j + 1 ] * data[ 2 * j + 1 ] );
      }
    output( bufferLength / 2 ) = 0.5 * fac * data[ 1 ] * data[ 1 ];

    for( int j = 0; j < bufferLength / 2; ++j )
      {
      matrix( j, i ) = output( j );
      }

    for( int j = bufferLength / 2; j < numberOfTimePoints; ++j )
      {
      matrix( j, i ) = 0;
      }
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( outputFileName.c_str() );
  writer->Update();

  return 0;
}
