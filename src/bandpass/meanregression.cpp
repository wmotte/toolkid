#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_inverse.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_vector.h"
#include "tkdCmdParser.h"

/** (c) Press et al., Numerical Recipes 3rd ed., p.784 **/
void fitline( const vnl_vector< float >& data, float& a, float& b, float& mean )
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
  mean = sy / static_cast< float >( numberOfPoints );
}

void removeTrend( vnl_vector< float >& signal, bool addMean = false )
{
  float a, b, mean;
  fitline( signal, a, b, mean );

  const int n = signal.size();
  for( int i = 0; i < n; ++i )
    {
    signal( i ) = signal( i ) - a * static_cast< float >( i ) - b;
    if ( addMean )
      {
      signal( i ) += mean;
      }
    }
}

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "meanregression", "Linear de-trending and global signal regression for 4-D fMRI datasets" );

  std::string inputFileName, outputFileName, maskFileName;
  bool detrend = false;
  bool addMean = false;

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

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetDescription( "3-D mask" );

  p.AddArgument( detrend, "detrend" )
    ->AddAlias( "d" )
    ->SetDescription( "Remove linear trend" );

  p.AddArgument( addMean, "add-mean" )
    ->AddAlias( "am" )
    ->SetDescription( "Add mean back to de-trended signal" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }


  typedef itk::Image< float, 4 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::Image< unsigned char, 3 > MaskType;
  typedef itk::ImageFileReader< MaskType > MaskReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  MaskReaderType::Pointer maskReader = MaskReaderType::New();

  reader->SetFileName( inputFileName.c_str() );
  reader->Update();

  ImageType::Pointer image = reader->GetOutput();
  reader = 0;

  float* buffer = image->GetPixelContainer()->GetBufferPointer();
  int timepoints = image->GetLargestPossibleRegion().GetSize()[ 3 ];
  int voxels = image->GetLargestPossibleRegion().GetNumberOfPixels() / timepoints;

  vnl_matrix_ref< float > matrix( timepoints, voxels, buffer );

  // construct mask
  vnl_vector< int > mask( voxels );
  mask.fill( 1 );

  if ( maskFileName != "" )
    {
    maskReader->SetFileName( maskFileName.c_str() );
    maskReader->Update();
    unsigned char* maskBuffer = maskReader->GetOutput()->GetPixelContainer()->GetBufferPointer();

    for( int i = 0; i < voxels; ++i )
      {
      mask[ i ] = ( maskBuffer[ i ] != 0 ) ? 1 : 0;
      }
    }

  // determine mean and remove linear trend
  vnl_vector< float > mean( timepoints );
  mean.fill( 0 );
  int count = 0;
  for( int i = 0; i < voxels; ++i )
    {
    vnl_vector< float > column = matrix.get_column( i );
    if ( detrend )
      {
      removeTrend( column, addMean );
      matrix.set_column( i, column );
      }

    if ( mask[ i ] == 0 )
      {
      continue;
      }

    mean += column;
    ++count;
    }

  // calculate mean signal
  mean /= static_cast< float >( count );

  // build masked input matrix
  vnl_matrix< float > B( timepoints, count );
  count = 0;
  for( int i = 0; i < voxels; ++i )
    {
    if ( mask[ i ] == 0 )
      {
      continue;
      }

    B.set_column( count++, matrix.get_column( i ) );
    }

  // regression (cf. Fox 2009, Appendix)
  vnl_matrix_ref< float > g( timepoints, 1, mean.data_block() );
//  vnl_matrix< float > gPlus = g.transpose() * ( 1.0 / dot_product( mean, mean ) );
//  vnl_matrix< float > BetaG = gPlus * B;
//  vnl_matrix< float > Bg = B - ( g * BetaG );

  /* GLM: Y = X Beta + epsilon
   *
   * Solve for Beta:
   *   Beta = ( X' X )^-1 X' Y
   *   pseudo-inverse of X'X with SVD,
   *     M  = U V  W', then
   *     M^-1 = W V^-1 U
   *     where V^-1 is reciprocal value for each diagonal element
   *
   * Then removing regressed variables:
   * Y' = Y - X Beta
   *    = Y - X ( X'X )^-1 X' Y
   */
  vnl_svd< float > svd( g.transpose() * g );
  vnl_matrix< float > Bg = B - g * svd.pinverse() * g.transpose() * B;

  // reconstruct
  count = 0;
  matrix.fill( 0 );
  for( int i = 0; i < voxels; ++i )
    {
    if ( mask[ i ] == 0 )
      {
      continue;
      }

    matrix.set_column( i, Bg.get_column( count++ ) );
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( image );
  writer->Update();

  return 0;
}

