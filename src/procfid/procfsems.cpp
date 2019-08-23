#include "itkImage.h"
#include <iostream>
#include "vnl/vnl_matrix.h"
#include "fidFID.h"
#include "fidFIDReader.h"
#include "fidFourier.h"
#include "tkdCmdParser.h"
#include "fidCommon.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "procfsems", "Fourier transform Varian FSEMS FID" );

  std::string inputFileName;
  std::string outputFileName;
  std::string outputRealFileName;
  std::string outputImaginaryFileName;
  bool outputKSpace = false;
  std::vector< int > zeroFill;

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
		->SetInput( "<int int>" )
		->SetDescription( "Zero-fill matrix" )
		->SetMinMax( 2, 2 );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef fid::FID::ComplexType ComplexType;
  typedef fid::FID::PrecisionType PrecisionType;
  typedef fid::Fourier< PrecisionType > FourierType;
  typedef vnl_matrix< ComplexType > MatrixType;
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
  int numberOfZ = fid->Procpar< int >( "ns" ); // slices
  int numberOfT = fid->GetNumberOfBlocks();
  int etl = fid->Procpar< int >( "etl" ); // echo train

  // zero-fill?
  int zeroFillRO = numberOfRO;
  int zeroFillPE = numberOfPE;

  if ( zeroFill.size() == 2 )
  	{
  	zeroFillRO = zeroFill[ 0 ];
  	zeroFillPE = zeroFill[ 1 ];
  	}

  vnl_matrix< int > table( numberOfPE / etl, etl );
  int* tableData = table.data_block();
  int tableMin = itk::NumericTraits< int >::max();

  if ( !fid->GetProcpar().Has( "pelist" ) )
    {
    std::string tablePath = fid::Common::GetAbsolutePath(
          argv[ 0 ],
          std::string( "tables/" ) + fid->GetProcpar().GetString( "petable" ) );

    std::ifstream in( tablePath.c_str() );
    std::string line;
    std::getline( in, line );
    for( int i = 0; i < numberOfPE / etl; ++i )
      {
      for( int j = 0; j < etl; ++j )
        {
        in >> table( i, j );
        tableMin = tableMin < table( i, j ) ? tableMin : table( i, j );
        }
      }
    }
  else
    {
    for( int i = 0; i < fid->GetProcpar().GetSize( "pelist" ); ++i )
      {
      tableData[ i ] = fid->Procpar< int >( "pelist", i );
      tableMin = tableMin < tableData[ i ] ? tableMin : tableData[ i ];
      }
    }

  // prepare output image

  OutputImageType::Pointer outputImage = OutputImageType::New();
  OutputImageType::SizeType size;
  OutputImageType::RegionType region;
  size[ 0 ] = zeroFillPE;
  size[ 1 ] = zeroFillRO;
  size[ 2 ] = numberOfZ;
  size[ 3 ] = numberOfT;
  region.SetSize( size );
  outputImage->SetRegions( region );
  outputImage->Allocate();

  ComplexType* outputPointer = outputImage->GetPixelContainer()->GetBufferPointer();

  // sort slices

  std::vector< float > pss = fid::Common::GetSlices( fid->GetProcpar() );
  std::vector< int > slices = fid::Common::GetSliceTable( pss );

  // calculate spacing and field of view

  float lro = fid->Procpar< float >( "lro" );
  float lpe = fid->Procpar< float >( "lpe" );

  float sum = 0;

  for( unsigned int i = 1; i < pss.size(); ++i )
    {
    sum += ( pss[ slices[ i ] ] - pss[ slices[ i - 1 ] ] );
    }

  float slcdist = pss.size() == 1 ? .1 : ( sum / static_cast< float >( pss.size() - 1 ) );

  OutputImageType::SpacingType spacing;
  spacing[ 0 ] = lpe / static_cast< float >( size[ 0 ] ) * 10.0;
  spacing[ 1 ] = lro / static_cast< float >( size[ 1 ] ) * 10.0;
  spacing[ 2 ] = slcdist * 10.0;
  spacing[ 3 ] = 1;

  outputImage->SetSpacing( spacing );

  float pss0 = fid->Procpar< float >( "pss0" ) * 10.0;

  OutputImageType::PointType origin;
  for( int i = 0; i < 3; ++i )
    {
    origin[ i ] = -0.5 * spacing[ i ] * static_cast< float >( size[ i ] );
    }
  origin[ 2 ] += pss0;
  origin[ 3 ] = 0;

  outputImage->SetOrigin( origin );

  // Fourier-transform

  int echoLength = numberOfPE / etl;

  FourierType fft( zeroFillPE, zeroFillRO );
  fft.SetOutputKSpace( outputKSpace );

  // reshape

  typedef itk::Image< ComplexType, 5 > ShapedDataType;
  ShapedDataType::Pointer shaped = ShapedDataType::New();
  ShapedDataType::RegionType shapedRegion;
  ShapedDataType::IndexType shapedIndex;
  ShapedDataType::SizeType shapedSize;

  shapedSize[ 0 ] = numberOfRO;
  shapedSize[ 1 ] = etl;
  shapedSize[ 2 ] = numberOfZ;
  shapedSize[ 3 ] = echoLength;
  shapedSize[ 4 ] = numberOfT;

  shapedIndex.Fill( 0 );

  shapedRegion.SetIndex( shapedIndex );
  shapedRegion.SetSize( shapedSize );

  shaped->SetRegions( shapedRegion );
  shaped->GetPixelContainer()->SetImportPointer(
  		fid->GetData()->GetPixelContainer()->GetBufferPointer(),
  		shapedRegion.GetNumberOfPixels(), false );

  shapedSize[ 1 ] = 1;
  shapedSize[ 2 ] = 1;
  shapedSize[ 3 ] = 1;
  shapedSize[ 4 ] = 1;

  // for each block
  for( int t = 0; t < numberOfT; ++t )
    {
    shapedIndex[ 4 ] = t;

    // for each slice
    for( int z = 0; z < numberOfZ; ++z )
      {
      // sorted slice
      shapedIndex[ 2 ] = slices[ z ];

      // fill slice matrix column-wise
      MatrixType sliceMatrix( numberOfPE, numberOfRO );

      for( int echoCount = 0; echoCount < echoLength; ++echoCount )
        {
        shapedIndex[ 3 ] = echoCount;

        for( int echo = 0; echo < etl; ++echo )
          {
          shapedIndex[ 1 ] = echo;

          // destination phase line
          int phase = table( echoCount, echo ) - tableMin;

          // copy

          ComplexType* pointer = &( shaped->GetPixel( shapedIndex ) );
          ComplexType* outputPointer = sliceMatrix[ phase ];
          for( int point = 0; point < numberOfRO; ++point, ++pointer, ++outputPointer )
          	{
          	( *outputPointer ) = ( *pointer );
          	}
          }
        }

      // Fourier-transform
      fft.ifft( sliceMatrix );

      for( unsigned int column = 0; column < sliceMatrix.cols(); ++column )
        {
        for( unsigned int row = 0; row < sliceMatrix.rows(); ++row, ++outputPointer )
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

