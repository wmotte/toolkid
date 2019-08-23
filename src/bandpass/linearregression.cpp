#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_inverse.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_vector.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "linearregression", "Linear regression with least-squares on 4-D fMRI datasets" );

  std::string inputFileName, outputFileName, maskFileName, betaFileName;
  std::vector< std::string > designFileNames;

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

  p.AddArgument( betaFileName, "output-beta" )
    ->AddAlias( "ob" )
    ->SetInput( "filename" )
    ->SetDescription( "Output 4D image with fit values" );

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetDescription( "3-D mask" );

  p.AddArgument( designFileNames, "design" )
    ->AddAlias( "d" )
    ->SetDescription( "Text files with design matrix (one column per file)" );

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

  int count = voxels;
  if ( maskFileName != "" )
    {
    maskReader->SetFileName( maskFileName.c_str() );
    maskReader->Update();
    unsigned char* maskBuffer = maskReader->GetOutput()->GetPixelContainer()->GetBufferPointer();

    count = 0;
    for( int i = 0; i < voxels; ++i )
      {
      mask[ i ] = ( maskBuffer[ i ] != 0 ) ? 1 : 0;
      count += mask[ i ];
      }
    }

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

  /* GLM: B = X Beta + epsilon
   *
   * Solve for Beta:
   *   Beta = ( X' X )^-1 X' B
   *   pseudo-inverse of X'X with SVD,
   *     M  = U V  W', then
   *     M^-1 = W V^-1 U
   *     where V^-1 is reciprocal value for each diagonal element
   *
   * Then removing regressed variables:
   * B' = B - X Beta
   *    = B - X ( X'X )^-1 X' B
   */

  int numberOfColumns = designFileNames.size();
  vnl_matrix< float > X( timepoints, numberOfColumns );

  for( int i = 0; i < numberOfColumns; ++i )
    {
    std::ifstream in( designFileNames[ i ].c_str() );
    for( int j = 0; j < timepoints; ++j )
      {
      in >> X( j, i );
      }
    }

  vnl_svd< float > svd( X.transpose() * X );
  vnl_matrix< float > Beta = svd.pinverse() * X.transpose() * B;
  vnl_matrix< float > fit = X * Beta;
  vnl_matrix< float > removed = B - fit;

  // prepare fitted output
  ImageType::Pointer outputFit = ImageType::New();
  outputFit->CopyInformation( image );

  ImageType::RegionType region = image->GetLargestPossibleRegion();
  outputFit->SetRegions( region );
  outputFit->Allocate();
  outputFit->FillBuffer( 0 );

  vnl_matrix_ref< float > fitMatrix( timepoints, voxels, outputFit->GetPixelContainer()->GetBufferPointer() );

  // prepare beta output
  ImageType::Pointer outputBeta = ImageType::New();
  outputBeta->CopyInformation( image );

  ImageType::SizeType size = region.GetSize();
  size[ 3 ] = numberOfColumns;
  region.SetSize( size );
  outputBeta->SetRegions( region );
  outputBeta->Allocate();
  outputBeta->FillBuffer( 0 );

  vnl_matrix_ref< float > betaMatrix( numberOfColumns, voxels, outputBeta->GetPixelContainer()->GetBufferPointer() );

  // reconstruct
  count = 0;
  matrix.fill( 0 );
  for( int i = 0; i < voxels; ++i )
    {
    if ( mask[ i ] == 0 )
      {
      continue;
      }

    matrix.set_column( i, removed.get_column( count ) );
    betaMatrix.set_column( i, Beta.get_column( count ) );
    fitMatrix.set_column( i, fit.get_column( count++ ) );
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( image );
  writer->Update();

  if ( betaFileName != "" )
    {
    writer->SetFileName( betaFileName.c_str() );
    writer->SetInput( outputBeta );
    writer->Update();
    }

  return 0;
}

