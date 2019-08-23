#include "dtifit.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "procparser.h"

namespace dtifit
{

struct DTIFit::Impl
{
  typedef itk::VectorImage< PixelType, 3 > ImageType;
  ImageType::Pointer Input;
  GradientTableType GradientTable;
  double BValue;
  OutputImageType::Pointer FA;
  OutputImageType::Pointer RA;
  OutputImageType::Pointer Trace;

  OutputImageType::Pointer L1;
  OutputImageType::Pointer L2;
  OutputImageType::Pointer L3;

  OutputVectorImageType::Pointer V1;
  OutputVectorImageType::Pointer V2;
  OutputVectorImageType::Pointer V3;

  TensorImageType::Pointer Tensor;

  RGBImageType::Pointer ColorFA;
};

DTIFit::DTIFit() : m_Impl( new Impl, true )
{
  m_Impl->BValue = 0;
}

DTIFit::~DTIFit()
{

}

void DTIFit::SetImage( InputImageType::Pointer input )
{
  // extract images to a vector image

  typedef itk::Image< PixelType, 3 > ImageType;
  typedef itk::ExtractImageFilter< InputImageType, ImageType > ExtractType;

  InputImageType::RegionType region = input->GetLargestPossibleRegion();
  InputImageType::IndexType index = region.GetIndex();
  InputImageType::SizeType size = region.GetSize();

  int numberOfImages = size[ 3 ];

  index.Fill( 0 );
  size[ 3 ] = 0;

  for( int i = 0; i < numberOfImages; ++i )
   {
   index[ 3 ] = i;
   region.SetIndex( index );
   region.SetSize( size );

   ExtractType::Pointer extract = ExtractType::New();
   extract->SetInput( input );
   extract->SetExtractionRegion( region );
   extract->Update();

   ImageType::Pointer output = extract->GetOutput();
   output->DisconnectPipeline();

   if ( !m_Impl->Input )
     {
     m_Impl->Input = Impl::ImageType::New();
     m_Impl->Input->CopyInformation( output );
     m_Impl->Input->SetRegions( output->GetLargestPossibleRegion() );
     m_Impl->Input->SetNumberOfComponentsPerPixel( numberOfImages );
     m_Impl->Input->Allocate();
     Impl::ImageType::PixelType vector( numberOfImages );
     m_Impl->Input->FillBuffer( vector );
     }

   itk::ImageRegionConstIterator< ImageType > itIn( output, output->GetLargestPossibleRegion() );
   itk::ImageRegionIterator< Impl::ImageType > itOut( m_Impl->Input, output->GetLargestPossibleRegion() );

   for( ; !itIn.IsAtEnd(); ++itIn, ++itOut )
     {
     Impl::ImageType::PixelType vector = itOut.Get();
     vector.SetElement( i, itIn.Value() );
     itOut.Set( vector );
     }
   }
}

void DTIFit::SetGradientTable( const GradientTableType& table )
{
  m_Impl->GradientTable = table;
}

void DTIFit::SetBValue( double bValue )
{
  m_Impl->BValue = bValue;
}

void DTIFit::Fit()
{
  // fit DTI

  typedef itk::DiffusionTensor3DReconstructionImageFilter< PixelType > FilterType;
  typedef FilterType::GradientDirectionContainerType ContainerType;
  typedef FilterType::GradientDirectionType GradientDirectionType;

  FilterType::Pointer filter = FilterType::New();
  ContainerType::Pointer table = ContainerType::New();

  // import gradient table
  table->Reserve( m_Impl->GradientTable.rows() );
  for( int i = 0; i < m_Impl->GradientTable.rows(); ++i )
    {
    GradientDirectionType direction = m_Impl->GradientTable.get_row( i );
    table->SetElement( i, direction );
    }

  filter->SetBValue( m_Impl->BValue );
  filter->SetGradientImage( table, m_Impl->Input );

  // set number of threads because of library bugs?
  filter->SetNumberOfThreads( 1 );
  filter->Update();

  // get output

  typedef FilterType::TensorImageType TensorImageType;
  typedef FilterType::TensorPixelType TensorPixelType;
  typedef TensorPixelType::EigenValuesArrayType EigenValuesArrayType;
  typedef TensorPixelType::EigenVectorsMatrixType EigenVectorsMatrixType;
  typedef TensorPixelType::RealValueType RealValueType;

  TensorImageType::Pointer output = filter->GetOutput();
  TensorImageType::RegionType region = output->GetLargestPossibleRegion();

  // compute eigen values

  OutputImageType::Pointer outputFA = OutputImageType::New();
  OutputImageType::Pointer outputRA = OutputImageType::New();
  OutputImageType::Pointer outputTrace = OutputImageType::New();

  OutputImageType::Pointer outputL1 = OutputImageType::New();
  OutputImageType::Pointer outputL2 = OutputImageType::New();
  OutputImageType::Pointer outputL3 = OutputImageType::New();

  RGBImageType::Pointer colorFA = RGBImageType::New();

  outputFA->CopyInformation( output );
  outputRA->CopyInformation( output );
  outputTrace->CopyInformation( output );

  outputL1->CopyInformation( output );
  outputL2->CopyInformation( output );
  outputL3->CopyInformation( output );

  colorFA->CopyInformation( output );

  outputFA->SetRegions( region );
  outputRA->SetRegions( region );
  outputTrace->SetRegions( region );

  outputL1->SetRegions( region );
  outputL2->SetRegions( region );
  outputL3->SetRegions( region );

  colorFA->SetRegions( region );

  outputFA->Allocate();
  outputRA->Allocate();
  outputTrace->Allocate();

  outputL1->Allocate();
  outputL2->Allocate();
  outputL3->Allocate();

  colorFA->Allocate();

  itk::ImageRegionConstIterator< TensorImageType > itTensor( output, region );

  itk::ImageRegionIterator< OutputImageType > itFA( outputFA, region );
  itk::ImageRegionIterator< OutputImageType > itRA( outputRA, region );
  itk::ImageRegionIterator< OutputImageType > itTrace( outputTrace, region );

  itk::ImageRegionIterator< OutputImageType > itL1( outputL1, region );
  itk::ImageRegionIterator< OutputImageType > itL2( outputL2, region );
  itk::ImageRegionIterator< OutputImageType > itL3( outputL3, region );

  itk::ImageRegionIterator< RGBImageType > itColorFA( colorFA, region );

  OutputVectorImageType::Pointer outputV1 = OutputVectorImageType::New();
  OutputVectorImageType::Pointer outputV2 = OutputVectorImageType::New();
  OutputVectorImageType::Pointer outputV3 = OutputVectorImageType::New();

  TensorImageType::SizeType tensorSize = region.GetSize();
  TensorImageType::SpacingType tensorSpacing = output->GetSpacing();
  TensorImageType::PointType tensorOrigin = output->GetOrigin();
  TensorImageType::DirectionType tensorDirection = output->GetDirection();

  OutputVectorImageType::RegionType vectorRegion;
  OutputVectorImageType::SizeType vectorSize;
  OutputVectorImageType::DirectionType vectorDirection;
  OutputVectorImageType::SpacingType vectorSpacing;
  OutputVectorImageType::PointType vectorOrigin;
  OutputVectorImageType::IndexType vectorIndex;

  vectorIndex.Fill( 0 );
  vectorDirection.SetIdentity();
  for( int i = 0; i < 3; ++i )
    {
    vectorSize[ i ] = tensorSize[ i ];
    vectorSpacing[ i ] = tensorSpacing[ i ];
    vectorOrigin[ i ] = tensorOrigin[ i ];

    for( int j = 0; j < 3; ++j )
      {
      vectorDirection[ i ][ j ] = tensorDirection[ i ][ j ];
      }
    }

  vectorSize[ 3 ] = 3;
  vectorSpacing[ 3 ] = 1;
  vectorOrigin[ 3 ] = 0;

  vectorRegion.SetIndex( vectorIndex );
  vectorRegion.SetSize( vectorSize );

  outputV1->SetRegions( vectorRegion );
  outputV1->SetSpacing( vectorSpacing );
  outputV1->SetOrigin( vectorOrigin );
  outputV1->SetDirection( vectorDirection );
  outputV1->Allocate();

  outputV2->CopyInformation( outputV1 );
  outputV2->SetRegions( vectorRegion );
  outputV2->Allocate();

  outputV3->CopyInformation( outputV1 );
  outputV3->SetRegions( vectorRegion );
  outputV3->Allocate();

  typedef itk::ImageRegionIterator< OutputVectorImageType > VectorIterator;
  std::vector< VectorIterator > vectorIterators;

  for( int i = 0; i < 3; ++i )
    {
    vectorSize[ 3 ] = 1;
    vectorIndex[ 3 ] = i;
    vectorRegion.SetSize( vectorSize );
    vectorRegion.SetIndex( vectorIndex );

    vectorIterators.push_back( VectorIterator( outputV3, vectorRegion ) );
    vectorIterators.push_back( VectorIterator( outputV2, vectorRegion ) );
    vectorIterators.push_back( VectorIterator( outputV1, vectorRegion ) );
    }

  for( ; !itTensor.IsAtEnd(); ++itTensor, ++itFA, ++itRA, ++itTrace, ++itL1, ++itL2, ++itL3, ++itColorFA )
    {
    const TensorPixelType& tensor = itTensor.Get();

    EigenValuesArrayType eigenValues;
    EigenVectorsMatrixType eigenVectors;
    tensor.ComputeEigenAnalysis( eigenValues, eigenVectors );

    itFA.Set( tensor.GetFractionalAnisotropy() );
    itRA.Set( tensor.GetRelativeAnisotropy() );
    itTrace.Set( tensor.GetTrace() );

    itL1.Set( eigenValues[ 2 ] );
    itL2.Set( eigenValues[ 1 ] );
    itL3.Set( eigenValues[ 0 ] );

    RGBPixelType& rgb = itColorFA.Value();
    rgb[ 0 ] = eigenValues[ 2 ] * 255. * itFA.Value();
    rgb[ 1 ] = eigenValues[ 1 ] * 255. * itFA.Value();
    rgb[ 2 ] = eigenValues[ 0 ] * 255. * itFA.Value();

    int iterator = 0;
    for( int i = 0; i < 3; ++i ) // for each component x, y, z
      {
      for( int j = 0; j < 3; ++j, ++iterator ) // for each
        {
        VectorIterator& it = vectorIterators[ iterator ];
        it.Set( eigenVectors( j, i ) );
        ++it;
        }
      }
    }

  m_Impl->FA = outputFA;
  m_Impl->RA = outputRA;
  m_Impl->Trace = outputTrace;

  m_Impl->L1 = outputL1;
  m_Impl->L2 = outputL2;
  m_Impl->L3 = outputL3;

  m_Impl->V1 = outputV1;
  m_Impl->V2 = outputV2;
  m_Impl->V3 = outputV3;

  m_Impl->Tensor = output;

  m_Impl->ColorFA = colorFA;
}

DTIFit::OutputImageType::Pointer DTIFit::GetFA()
{
  return m_Impl->FA;
}

DTIFit::OutputImageType::Pointer DTIFit::GetRA()
{
  return m_Impl->RA;
}

DTIFit::OutputImageType::Pointer DTIFit::GetTrace()
{
  return m_Impl->Trace;
}

DTIFit::OutputImageType::Pointer DTIFit::GetL1()
{
  return m_Impl->L1;
}

DTIFit::OutputImageType::Pointer DTIFit::GetL2()
{
  return m_Impl->L2;
}

DTIFit::OutputImageType::Pointer DTIFit::GetL3()
{
  return m_Impl->L3;
}

DTIFit::OutputVectorImageType::Pointer DTIFit::GetV1()
{
  return m_Impl->V1;
}

DTIFit::OutputVectorImageType::Pointer DTIFit::GetV2()
{
  return m_Impl->V2;
}

DTIFit::OutputVectorImageType::Pointer DTIFit::GetV3()
{
  return m_Impl->V3;
}

DTIFit::TensorImageType::Pointer DTIFit::GetTensor()
{
  return m_Impl->Tensor;
}

DTIFit::RGBImageType::Pointer DTIFit::GetColorFA()
{
  return m_Impl->ColorFA;
}

void DTIFit::Run( InputImageType::Pointer input, const Procparser& pp, const std::string& outputFileName, const std::string& outputExtension, bool mirrorGradientTable )
{
  int numberOfRO = pp.GetSize( "dro" );
  int numberOfImages = input->GetLargestPossibleRegion().GetSize()[ 3 ];

  vnl_matrix< double > table( numberOfRO, 3 );

  for( int i = 0; i < numberOfRO; ++i )
    {
    table( i, 0 ) = pp.GetAs< double >( "dpe", i );
    table( i, 1 ) = pp.GetAs< double >( "dro", i );
    table( i, 2 ) = pp.GetAs< double >( "dsl", i );
    }

  if ( numberOfImages == 2 * numberOfRO )
    {
    vnl_matrix< double > twice( numberOfRO * 2, 3 );

    if ( !mirrorGradientTable )
      {
      for( int i = 0; i < numberOfRO; ++i )
        {
        twice.set_row( i * 2, table.get_row( i ) );
        twice.set_row( i * 2 + 1, table.get_row( i ) * -1. );
        }
      }
    else
      {
      for( int i = 0; i < numberOfRO; ++i )
        {
        twice.set_row( i, table.get_row( i ) );
        }

      for( int i = 0; i < numberOfRO; ++i )
        {
        twice.set_row( i + numberOfRO, table.get_row( i ) * -1. );
        }
      }

    table = twice;
    numberOfRO = 2 * numberOfRO;
    }

/*
  if ( numberOfImages != numberOfRO )
    {
    std::cout << "WARNING: number of acquired images (" << numberOfImages << ") != number of gradient directions (" << numberOfRO << ")" << std::endl;

    if ( numberOfImages > numberOfRO )
      {
      return;
      }

    table = table.extract( numberOfImages, 3 );
    numberOfRO = numberOfImages;
    }
*/
  std::string orient = pp.GetAs< std::string >( "orient" );

  if ( orient == "cor" )
    {
    vnl_matrix< double > permute( 3, 3 );
    permute.set_identity();

    permute( 0, 0 ) = -1;
    permute( 0, 1 ) = 0;
    permute( 0, 2 ) = 0;

    permute( 1, 0 ) = 0;
    permute( 1, 1 ) = 0;
    permute( 1, 2 ) = 1;

    permute( 2, 0 ) = 0;
    permute( 2, 1 ) = 1;
    permute( 2, 2 ) = 0;

    table = ( permute * table.transpose() ).transpose();
    }
  else if ( orient == "trans90" )
    {
    vnl_matrix< double > permute( 3, 3 );
    permute.set_identity();

    permute( 0, 0 ) = 0;
    permute( 0, 1 ) = 1;
    permute( 0, 2 ) = 0;

    permute( 1, 0 ) = -1;
    permute( 1, 1 ) = 0;
    permute( 1, 2 ) = 0;

    permute( 2, 0 ) = 0;
    permute( 2, 1 ) = 0;
    permute( 2, 2 ) = 1;

    table = ( permute * table.transpose() ).transpose();
    }
  else if ( orient == "trans" )
    {
    vnl_matrix< double > permute( 3, 3 );
    permute.set_identity();

    permute( 0, 0 ) = 1;
    permute( 0, 1 ) = 0;
    permute( 0, 2 ) = 0;

    permute( 1, 0 ) = 0;
    permute( 1, 1 ) = -1;
    permute( 1, 2 ) = 0;

    permute( 2, 0 ) = 0;
    permute( 2, 1 ) = 0;
    permute( 2, 2 ) = 1;

    table = ( permute * table.transpose() ).transpose();
    }
  else
    {
    std::cout << "WARNING: unsupported orientation: " << orient << std::endl;
    }

  double gdiff = pp.GetAs< double >( "gdiff" );
  double tdelta = pp.GetAs< double >( "tdelta" );
  double tDELTA = pp.GetAs< double >( "tDELTA" );

  //Price and Kuchel, J Magnetic Resonance 94 (1991) 133-139 eq.5 page 137!
  double bValue =
        ( 2. / vnl_math::pi ) * ( 2. / vnl_math::pi )
      * ( tdelta * tdelta )
      * ( gdiff * gdiff )
      * ( tDELTA - tdelta / 4. ) // 4 not 3 !!
      * ( 26751.98775 * 26751.98775 * 0.01 );

  this->SetBValue( bValue );
  this->SetImage( input );
  this->SetGradientTable( table );

  std::ofstream bvals( ( outputFileName + std::string( "_bvals" ) ).c_str() );
  std::ofstream bvecs( ( outputFileName + std::string( "_bvecs" ) ).c_str() );

  for( int i = 0; i < numberOfRO; ++i )
    {
    vnl_vector< double > vector = table.get_row( i );
    bvals << ( vector( 0 ) * vector( 0 ) + vector( 1 ) * vector( 1 ) + vector( 2 ) * vector( 2 ) ) * bValue << std::endl;
    bvecs << vector( 0 ) << " " << vector( 1 ) << " " << vector( 2 ) << std::endl;
    }

  this->Fit();

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typedef itk::ImageFileWriter< OutputVectorImageType > VectorWriterType;
  typedef itk::ImageFileWriter< TensorImageType > TensorWriterType;
  typedef itk::ImageFileWriter< RGBImageType > RGBWriterType;

  WriterType::Pointer writer = WriterType::New();
  VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
  TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
  RGBWriterType::Pointer rgbWriter = RGBWriterType::New();

  writer->SetFileName( ( outputFileName + std::string( "_FA." ) + outputExtension ).c_str() );
  writer->SetInput( this->GetFA() );
  writer->Update();

  writer->SetFileName( ( outputFileName + std::string( "_RA." ) + outputExtension ).c_str() );
  writer->SetInput( this->GetRA() );
  writer->Update();

  writer->SetFileName( ( outputFileName + std::string( "_Trace." ) + outputExtension ).c_str() );
  writer->SetInput( this->GetTrace() );
  writer->Update();

  writer->SetFileName( ( outputFileName + std::string( "_L1." ) + outputExtension ).c_str() );
  writer->SetInput( this->GetL1() );
  writer->Update();

  writer->SetFileName( ( outputFileName + std::string( "_L2." ) + outputExtension ).c_str() );
  writer->SetInput( this->GetL2() );
  writer->Update();

  writer->SetFileName( ( outputFileName + std::string( "_L3." ) + outputExtension ).c_str() );
  writer->SetInput( this->GetL3() );
  writer->Update();

  vectorWriter->SetFileName( ( outputFileName + std::string( "_V1." ) + outputExtension ).c_str() );
  vectorWriter->SetInput( this->GetV1() );
  vectorWriter->Update();

  vectorWriter->SetFileName( ( outputFileName + std::string( "_V2." ) + outputExtension ).c_str() );
  vectorWriter->SetInput( this->GetV2() );
  vectorWriter->Update();

  vectorWriter->SetFileName( ( outputFileName + std::string( "_V3." ) + outputExtension ).c_str() );
  vectorWriter->SetInput( this->GetV3() );
  vectorWriter->Update();

  tensorWriter->SetFileName( ( outputFileName + std::string( "_Tensor." ) + outputExtension ).c_str() );
  tensorWriter->SetInput( this->GetTensor() );
  tensorWriter->Update();

  rgbWriter->SetFileName( ( outputFileName + std::string( "_ColorFA." ) + outputExtension ).c_str() );
  rgbWriter->SetInput( this->GetColorFA() );
  rgbWriter->Update();
}

} // end namespace dtifit
