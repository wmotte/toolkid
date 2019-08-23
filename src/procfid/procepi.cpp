#include "itkImage.h"
#include <iostream>
#include <fstream>
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"
#include "fidFID.h"
#include "fidFIDReader.h"
#include "fidFourier.h"
#include "itkImageRegionIterator.h"
#include "tkdCmdParser.h"
#include "fidEPITools.h"
#include "fidCommon.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "procepi", "Fourier transform single-shot EPI FID" );

  std::string inputFileName;
  std::string outputFileName;
  std::string referenceFileName;
  std::string outputRealFileName;
  std::string outputImaginaryFileName;
  bool outputKSpace = false;
  std::vector< int > zeroFill;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "FID input folder" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
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

  p.AddArgument( referenceFileName, "reference" )
    ->AddAlias( "r" )
    ->SetDescription( "FID reference folder" )
    ->SetRequired( true );

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
  typedef fid::FID::BlockType BlockType;
  typedef vnl_matrix< ComplexType > MatrixType;
  typedef itk::Image< ComplexType, 4 > OutputImageType;

  // reference image

  fid::FIDReader::Pointer referenceReader = fid::FIDReader::New();
  referenceReader->SetFileName( referenceFileName );

  try
    {
    referenceReader->Read();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Error reading reference FID: " << e.GetDescription() << std::endl;
    return -1;
    }

  fid::FID::Pointer reference = referenceReader->GetFID();
  std::vector< int > referenceSlices =
    fid::Common::GetSliceTable(
        fid::Common::GetSlices( reference->GetProcpar() ) );

  // load phase-encoding table

  vnl_matrix< int > table;
  std::string tablePath = fid::Common::GetAbsolutePath(
      argv[ 0 ],
      std::string( "tables/" ) + reference->Procpar< std::string >( "ropat" ) + ".txt" );

  if ( !fid::Common::LoadTable( tablePath, table ) )
    {
    std::cerr << "Can not open table: " << tablePath << std::endl;
    return -1;
    }

  int np = table.rows();
  int nv = table.cols();

  int numberOfRO = pow( 2, ceil( log( np ) / log( 2. ) ) );
  int numberOfPE = pow( 2, ceil( log( nv ) / log( 2. ) ) );
  int numberOfZ = reference->Procpar< int >( "ns" ); // slices

  // zero-fill?
  int zeroFillRO = numberOfRO;
  int zeroFillPE = numberOfPE;

  if ( zeroFill.size() == 2 )
  	{
  	zeroFillRO = zeroFill[ 0 ];
  	zeroFillPE = zeroFill[ 1 ];
  	}

  // rearrange reference FID

  vnl_matrix< ComplexType > referenceMatrix( numberOfZ, nv * np );

  for( int i = 0; i < reference->GetNumberOfBlocks(); ++i )
    {
    BlockType block = reference->GetBlock( i );

    for( int j = 0; j < nv; ++j )
      {
      for( int k = 0; k < np; ++k )
        {
				referenceMatrix.set_column( j * np + k, block.get_column( table( k, j ) ) );
				}
			}
    }

  // build phasemaps (for every slice: numberOfPE x numberOfRO)

  std::vector< vnl_matrix< PrecisionType > > phaseMaps;

  for( int z = 0; z < numberOfZ; ++z )
    {
    vnl_matrix< PrecisionType > phaseMap( numberOfPE, zeroFillRO );
    phaseMap.fill( itk::NumericTraits< PrecisionType >::Zero );

    // get slice
    vnl_vector< ComplexType > row = referenceMatrix.get_row( z );

    // import as matrix
    vnl_matrix_ref< ComplexType > sliceMatrix1( nv, np, row.data_block() );
    vnl_matrix< ComplexType > sliceMatrix( zeroFillRO, numberOfPE );

    for( int y = 0; y < numberOfPE; ++y )
      {
      fid::Fourier< PrecisionType > fft( zeroFillRO );
      vnl_vector< ComplexType > row = sliceMatrix1.get_row( y );
      fft.ifft( row );
      sliceMatrix.set_column( y, row );
      }

    // compute weights for fit

    vnl_vector< PrecisionType > weights = fid::EPITools::CalculateWeights( sliceMatrix );

    for( int y = 0; y < numberOfPE; ++y )
      {
      phaseMap.set_row( y, fid::EPITools::BuildPhaseMap( sliceMatrix.get_column( y ), weights ) );
      }

    phaseMaps.push_back( phaseMap );
    }

  // read image FID

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

  // for Annette's single-shot EPI:
  // every block is an acquired volume
  // every trace is a slice
  // every point needs to be reordered according to a read-out pattern matrix

  int numberOfT = fid->GetNumberOfBlocks();

  // prepare output image

  OutputImageType::Pointer outputImage = OutputImageType::New();
  OutputImageType::RegionType region;

  OutputImageType::SizeType size;
  size[ 0 ] = zeroFillPE;
  size[ 1 ] = zeroFillRO;
  size[ 2 ] = numberOfZ;
  size[ 3 ] = numberOfT;

  OutputImageType::IndexType index;
  index[ 0 ] = 0;
  index[ 1 ] = 0;
  index[ 2 ] = 0;
  index[ 3 ] = 0;

  region.SetIndex( index );
  region.SetSize( size );
  outputImage->SetRegions( region );
  outputImage->Allocate();

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

  // for each block

  for( int t = 0; t < numberOfT; ++t )
    {
    BlockType block = fid->GetBlock( t );
    vnl_matrix< ComplexType > fidMatrix( block.rows(), nv * np );
    fidMatrix.fill( itk::NumericTraits< ComplexType >::Zero );

    for( int j = 0; j < nv; ++j )
      {
      for( int k = 0; k < np; ++k )
        {
        fidMatrix.set_column( j * np + k, block.get_column( table( k, j ) ) );
        }
      }

    // for each slice
    for( int z = 0; z < numberOfZ; ++z)
      {
      vnl_matrix_ref< ComplexType > sliceMatrix( nv, np, fidMatrix.data_array()[ slices[ z ] ] );

      // get slice's phasemap
      const vnl_matrix< PrecisionType >& phaseMap = phaseMaps[ referenceSlices[ z ] ];

      // intermediate matrix
      vnl_matrix< ComplexType > intermediate( zeroFillRO, numberOfPE );

      for( int y = 0; y < numberOfPE; ++y )
        {
        vnl_vector< ComplexType > trace = sliceMatrix.get_row( y );
        fid::Fourier< PrecisionType > fft1( zeroFillRO );
//        fft1.SetOutputKSpace( outputKSpace );
        fft1.ifft( trace );

        // phase map
        fid::EPITools::PhaseWithMap( trace, phaseMap.get_row( y ) );

        if ( outputKSpace )
          {
          fft1.ifft( trace );
          }

        // store
        intermediate.set_column( y, trace );
        }

      vnl_matrix< ComplexType > intermediate2( zeroFillRO, zeroFillPE );

      for( int x = 0; x < zeroFillRO; ++x )
        {
        vnl_vector< ComplexType > trace = intermediate.get_row( x );
        fid::Fourier< PrecisionType > fft1( zeroFillPE );
        fft1.SetOutputKSpace( outputKSpace );
        fft1.ifft( trace );

        intermediate2.set_row( x, trace );
        }

      // save to output
      index[ 2 ] = z;
      index[ 3 ] = t;

      size[ 2 ] = 1;
      size[ 3 ] = 1;

      region.SetIndex( index );
      region.SetSize( size );

      itk::ImageRegionIterator< OutputImageType > it( outputImage, region );

      for( int row = 0, column = 0; !it.IsAtEnd(); ++it, ++column )
        {
        if ( column >= zeroFillPE )
          {
          column = 0;
          ++row;
          }

        it.Set( intermediate2( row, column ) );
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
