#include "itkImage.h"
#include "vnl/vnl_matrix.h"
#include "fidFIDReader.h"
#include "fidFID.h"
#include "fidFourier.h"
#include "fidCommon.h"
#include "itkImageRegionIterator.h"
#include "tkdCmdParser.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "mrfitCommon.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "procmsme", "Fourier transform Varian MSME FID" );

  std::string inputFileName;
  std::string outputFileName;
  std::vector< int > zeroFill;
  std::string outputRealFileName;
  std::string outputImaginaryFileName;
  std::string outputFileNameT2;
  int algorithm = 0;
  double maxT2 = 10;
  bool outputKSpace = false;

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

  p.AddArgument( zeroFill, "zerofill" )
		->AddAlias( "zf" )
		->SetInput( "<int int>" )
		->SetDescription( "Zero-fill matrix" )
		->SetMinMax( 2, 2 );

  p.AddArgument( outputFileNameT2, "output-T2map" )
    ->AddAlias( "t2map" )
    ->SetDescription( "Perform T2-mapping" );

  p.AddArgument( algorithm, "algorithm" )
    ->AddAlias( "a" )
    ->SetDescription( "Fitting algorithm: 0=LINEAR, 1=NON_LINEAR, 2=NON_LINEAR_WITH_CONSTANT (default: 0)" );

  p.AddArgument( maxT2, "maximum-T2" )
    ->AddAlias( "max" )
    ->SetDescription( "Maximum T2 (default: 10)" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef fid::FID::PrecisionType PrecisionType;
  typedef fid::FID::ComplexType ComplexType;
  typedef fid::FID::DataType DataType;
  typedef fid::Fourier< PrecisionType > FourierType;
  typedef vnl_matrix< ComplexType > MatrixType;
  typedef itk::Image< ComplexType, 4 > OutputImageType;

  // Read FID and procpar

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

  // permute FID

  typedef itk::PermuteAxesImageFilter< DataType > PermuteDataType;
  PermuteDataType::Pointer permuteData = PermuteDataType::New();
  permuteData->SetInput( fid->GetData() );
  PermuteDataType::PermuteOrderArrayType dataPermutation;
  dataPermutation[ 0 ] = 0;
  dataPermutation[ 1 ] = 2;
  dataPermutation[ 2 ] = 1;
  permuteData->SetOrder( dataPermutation );
  permuteData->Update();

  // determine dimensions
  int numberOfRO = fid->GetNumberOfPoints();
  int numberOfPE = fid->Procpar< int >( "nv" );
  int numberOfEchos = fid->Procpar< int >( "ne" );
  int numberOfZ = fid->Procpar< int >( "ns" );

  // reshape
  typedef itk::Image< ComplexType, 4 > ProcessImageType;
  ProcessImageType::Pointer processImage = ProcessImageType::New();
  ProcessImageType::RegionType processRegion;
  ProcessImageType::SizeType processSize;
  processSize[ 0 ] = numberOfRO;
  processSize[ 1 ] = numberOfPE;
  processSize[ 2 ] = numberOfEchos;
  processSize[ 3 ] = numberOfZ;
  processRegion.SetSize( processSize );
  processImage->SetRegions( processRegion );
  processImage->GetPixelContainer()->SetImportPointer(
  		permuteData->GetOutput()->GetPixelContainer()->GetBufferPointer(),
  		processRegion.GetNumberOfPixels(), false );

  // zero-fill?
  int zeroFillRO = numberOfRO;
  int zeroFillPE = numberOfPE;

  if ( zeroFill.size() == 2 )
  	{
  	zeroFillRO = zeroFill[ 0 ];
  	zeroFillPE = zeroFill[ 1 ];
  	}

  // prepare output image

  OutputImageType::Pointer outputImage = OutputImageType::New();
  OutputImageType::RegionType region;

  OutputImageType::SizeType size;
  size[ 0 ] = zeroFillPE;
  size[ 1 ] = zeroFillRO;
  size[ 2 ] = numberOfZ;
  size[ 3 ] = numberOfEchos;

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
  spacing[ 0 ] = lpe / static_cast< float >( zeroFillPE ) * 10.0;
  spacing[ 1 ] = lro / static_cast< float >( zeroFillRO ) * 10.0;
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
  ProcessImageType::IndexType processIndex;
  processIndex.Fill( 0 );

  processSize[ 2 ] = 1;
  processSize[ 3 ] = 1;

  // for each block
  for( int t = 0; t < numberOfEchos; ++t )
    {
    processIndex[ 2 ] = t;

    // for each slice
    for( int z = 0; z < numberOfZ; ++z )
      {
      // sorted slice
      int currentSlice = slices[ z ];
      processIndex[ 3 ] = currentSlice;

      // fill slice matrix column-wise
      MatrixType sliceMatrix( numberOfPE, numberOfRO );
      FourierType        fft( zeroFillPE, zeroFillRO );

      fft.SetOutputKSpace( outputKSpace );
      // for each row in FID, get trace corresponding to current slice
      processRegion.SetIndex( processIndex );
      processRegion.SetSize( processSize );

      // copy

      itk::ImageRegionConstIterator< ProcessImageType > processIt( processImage, processRegion );
      ComplexType* sliceMatrixPointer = sliceMatrix.data_block();

      for( int sliceMatrixSize = processRegion.GetNumberOfPixels();
					 sliceMatrixSize > 0;
					 ++sliceMatrixPointer, ++processIt, --sliceMatrixSize )
        {
        ( *sliceMatrixPointer ) = processIt.Value();
        }

      // flip

      if ( t % 2 != 0 )
      	{
      	sliceMatrix.flipud();
      	}

      // Fourier-transform
      fft.ifft( sliceMatrix );

      // save to output

      index[ 2 ] = z;
      index[ 3 ] = t;

      size[ 2 ] = 1;
      size[ 3 ] = 1;

      region.SetIndex( index );
      region.SetSize( size );

      itk::ImageRegionIterator< OutputImageType > it( outputImage, region );
      for( unsigned int column = 0; column < sliceMatrix.cols(); ++column )
        {
        for( unsigned int row = 0; row < sliceMatrix.rows(); ++row, ++it )
          {
          it.Set( sliceMatrix( row, column ) );
          }
        }
      }
    }

  // reorient data
  outputImage = fid::Common::ReorientCoronal( outputImage, fid->GetProcpar() );

  // write results

  itk::Image< float, 4 >::Pointer magnitude = fid::Common::ComplexToMagnitude( outputImage );
  fid::Common::SaveSeries( magnitude, outputFileName );

  if ( outputRealFileName != "" )
    {
    fid::Common::SaveSeries(
        fid::Common::ComplexToReal( outputImage ),
        outputRealFileName );
    }

  if ( outputImaginaryFileName != "" )
    {
    fid::Common::SaveSeries(
        fid::Common::ComplexToImaginary( outputImage ),
        outputImaginaryFileName );
    }

  // t2 mapping
  if ( outputFileNameT2 != "" )
    {
    std::vector< double > echos;
    if ( echos.size() == 0 )
      {
      const Procparser& pp = fid->GetProcpar();
      float te = fid->Procpar< float >( "te" );
      int ne = fid->Procpar< int >( "ne" );

      for( int i = 0; i < ne; ++i )
        {
        echos.push_back( te * static_cast< float >( i +1 ) );
        }
      }

    mrfit::MRFit::Pointer fit = mrfit::MRFit::New();
    fit->SetInput( magnitude );
    fit->SetTimes( echos );
    fit->FitT2( static_cast< mrfit::MRFit::T2FittingType >( algorithm ), maxT2 );

    fid::Common::SaveImage( fit->GetMap( 0 ), outputFileNameT2 );
    }


  return 0;
}
