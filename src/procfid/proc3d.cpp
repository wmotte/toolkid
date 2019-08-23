#include "itkImage.h"
#include "fidFID.h"
#include "fidFIDReader.h"
#include "fidFourier.h"
#include "tkdCmdParser.h"
#include "fidCommon.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "proc3d", "Fourier transform Varian FSE3D/GE3D FID" );

  std::string inputFileName;
  std::string outputFileName;
  std::string outputRealFileName;
  std::string outputImaginaryFileName;
  std::vector< int > zeroFill;
  bool outputKSpace = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "path" )
    ->SetDescription( "FID input folder" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "path" )
    ->SetDescription( "Output image" )
    ->SetRequired( true );

  p.AddArgument( outputRealFileName, "output-real" )
    ->AddAlias( "or" )
    ->SetDescription( "Output real image" );

  p.AddArgument( outputImaginaryFileName, "output-imag" )
    ->AddAlias( "oi" )
    ->SetDescription( "Output imaginary image" );

  p.AddArgument( outputKSpace, "k-space" )
    ->AddAlias( "k" )
    ->SetDescription( "Output k-space (no Fourier transform)" );

  p.AddArgument( zeroFill, "zerofill" )
		->AddAlias( "zf" )
		->SetInput( "<int int int>" )
		->SetDescription( "Zero-fill matrix" )
		->SetMinMax( 3, 3 );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef fid::FID::PrecisionType PrecisionType;
  typedef fid::FID::ComplexType ComplexType;
  typedef vnl_matrix< ComplexType > MatrixType;
  typedef fid::Fourier< PrecisionType > FourierType;
  typedef itk::Image< ComplexType, 4 > OutputImageType;

  fid::FIDReader::Pointer reader = fid::FIDReader::New();
  reader->SetFileName( inputFileName );

  try
    {
    reader->Read();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Error reading FID: " << e.GetDescription() << std::endl;
    return -1;
    }

  fid::FID::Pointer fid = reader->GetFID();

  int numberOfRO = fid->GetNumberOfPoints();
  int numberOfPE = fid->Procpar< int >( "nv" ); // columns
  int numberOfPE2 = fid->Procpar< int >( "nv2" ); // slices
  int numberOfBlocks = fid->GetNumberOfBlocks();

  // zero-fill?
  int zeroFillRO = numberOfRO;
  int zeroFillPE = numberOfPE;
  int zeroFillPE2 = numberOfPE2;

  if ( zeroFill.size() == 3 )
  	{
  	zeroFillRO = zeroFill[ 0 ];
  	zeroFillPE = zeroFill[ 1 ];
  	zeroFillPE2 = zeroFill[ 2 ];
  	}

  // prepare output image

  OutputImageType::Pointer outputImage = OutputImageType::New();
  OutputImageType::SizeType size;
  OutputImageType::RegionType region;
  size[ 0 ] = zeroFillPE;
  size[ 1 ] = zeroFillRO;
  size[ 2 ] = zeroFillPE2;
  size[ 3 ] = numberOfBlocks;
  region.SetSize( size );
  outputImage->SetRegions( region );
  outputImage->Allocate();

  ComplexType* outputPointer = outputImage->GetPixelContainer()->GetBufferPointer();

  // calculate spacing and field of view

  float lro  = fid->Procpar< float >( "lro" );
  float lpe  = fid->Procpar< float >( "lpe" );
  float lpe2 = fid->Procpar< float >( "lpe2" );

  OutputImageType::SpacingType spacing;
  spacing[ 0 ] = lpe  / static_cast< float >( size[ 0 ] ) * 10.0;
  spacing[ 1 ] = lro  / static_cast< float >( size[ 1 ] ) * 10.0;
  spacing[ 2 ] = lpe2 / static_cast< float >( size[ 2 ] ) * 10.0;
  spacing[ 3 ] = 1;

  outputImage->SetSpacing( spacing );

  OutputImageType::PointType origin;
  for( int i = 0; i < 3; ++i )
    {
    origin[ i ] = -0.5 * spacing[ i ] * static_cast< float >( size[ i ] );
    }
  origin[ 3 ] = 0;

  outputImage->SetOrigin( origin );

  // Fourier-transform

  // FID: trace x phase x slice x array
  FourierType fft1( zeroFillPE2 );
  FourierType fft2( zeroFillPE, zeroFillRO );

  fft1.SetOutputKSpace( outputKSpace );
  fft2.SetOutputKSpace( outputKSpace );

  // for each block
  for( int t = 0; t < numberOfBlocks; ++t )
    {
    fid::FID::BlockType block = fid->GetBlock( t );
    MatrixType slices( numberOfPE2, zeroFillRO * zeroFillPE );
    slices.fill( itk::NumericTraits< ComplexType >::Zero );

    // for each slice
    for( int z = 0; z < numberOfPE2; ++z )
      {
      // fill slice matrix column-wise
      MatrixType sliceMatrix( numberOfPE, numberOfRO );

      // for each row in FID, get trace corresponding to current slice
      for( int y = 0; y < numberOfPE; ++y )
        {
        sliceMatrix.set_row( y, block.get_row( numberOfPE * z + y ) );
        }

      fft2.ifft( sliceMatrix );
      slices.set_row( z, sliceMatrix.data_block() );
      }

    // 1-D Fourier transform
    MatrixType slices2( zeroFillPE2, zeroFillRO * zeroFillPE );
    for( int i = 0; i < ( zeroFillRO * zeroFillPE ); ++i )
      {
      vnl_vector< ComplexType > sliceVector = slices.get_column( i );
      fft1.ifft( sliceVector );
      slices2.set_column( i, sliceVector );
      }
    slices = slices2;

    // fill output (shift center)
    for( int i = zeroFillPE2 - 1; i >= 0; --i )
      {
      vnl_vector< ComplexType > sliceVector = slices.get_row( i );
      vnl_matrix_ref< ComplexType > sliceMatrix( zeroFillPE, zeroFillRO, sliceVector.data_block() );

      for( int column = 0; column < zeroFillRO; ++column )
        {
        for( int row = 0; row < zeroFillPE; ++row, ++outputPointer )
          {
          ( *outputPointer ) = sliceMatrix( row, column );
          }
        }
      }
    }

  // reorient data
  outputImage = fid::Common::ReorientCoronal( outputImage, fid->GetProcpar() );

  // write results
  if ( outputFileName != "" )
    {
    fid::Common::SaveSeries(
        fid::Common::ComplexToMagnitude( outputImage ), outputFileName );
    }

  if ( outputRealFileName != "" )
    {
    fid::Common::SaveSeries(
        fid::Common::ComplexToReal( outputImage ), outputRealFileName );
    }

  if ( outputImaginaryFileName != "" )
    {
    fid::Common::SaveSeries(
        fid::Common::ComplexToImaginary( outputImage ), outputImaginaryFileName );
    }

  return 0;
}

