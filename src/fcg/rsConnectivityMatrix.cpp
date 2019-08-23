#include "rsConnectivityMatrix.h"
#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"
#include "rsCommon.h"
#include "rsPCA.h"

namespace rs
{

ConnectivityMatrix::ConnectivityMatrix()
: m_SkipInitial( 0 ), m_SkipEnd( 0 )
{
  m_LabelMin = 0;
  m_LabelMax = 0;
}

ConnectivityMatrix::~ConnectivityMatrix()
{
}

void ConnectivityMatrix::SetInput( InputType::Pointer input )
{
  m_Input = input;
}

void ConnectivityMatrix::SetLabels( LabelType::Pointer labels )
{
  m_Labels = labels;
}

void ConnectivityMatrix::SetSkipTimePoints( int initial, int end )
{
  m_SkipInitial = initial;
  m_SkipEnd = end;
}

void ConnectivityMatrix::SetLabelRange( LabelPixelType min, LabelPixelType max )
{
  m_LabelMin = min;
  m_LabelMax = max;
}

void ConnectivityMatrix::Initialize()
{
  // get number of voxels and time points
  InputType::RegionType inputRegion = m_Input->GetLargestPossibleRegion();
  InputType::SizeType inputSize = inputRegion.GetSize();
  m_NumberOfVoxels = inputSize[ 0 ] * inputSize[ 1 ] * inputSize[ 2 ];
  m_NumberOfSamples = inputSize[ 3 ];

  // determine number of labels
  m_LabelSet.clear();
  int labelCounter = 0;

  if ( m_LabelMin != 0 && m_LabelMax != 0 && m_LabelMax >= m_LabelMin )
    {
    for( LabelPixelType i = m_LabelMin; i <= m_LabelMax; ++i )
      {
      m_LabelSet[ i ] = labelCounter++;
      }
    }
  else
    {
    itk::ImageRegionConstIterator< LabelType > it( m_Labels, m_Labels->GetLargestPossibleRegion() );
  //  int labelMax = 1;
    for( ; !it.IsAtEnd(); ++it )
      {
      if ( it.Value() < 1 )
        {
        continue;
        }

      if ( m_LabelSet.find( it.Value() ) == m_LabelSet.end() )
        {
        m_LabelSet[ it.Value() ] = labelCounter++;
        }

  //    labelMax = it.Value() > labelMax ? it.Value() : labelMax;
      }
    }

  m_NumberOfLabels = labelCounter;

  // initialize output matrix
  m_Output = OutputType::New();
  OutputType::RegionType outputRegion;
  OutputType::SizeType size;
  size[ 0 ] = m_NumberOfLabels;
  size[ 1 ] = size[ 0 ];
  outputRegion.SetSize( size );
  m_Output->SetRegions( outputRegion );
  m_Output->Allocate();
  m_Output->FillBuffer( itk::NumericTraits< PixelType >::Zero );
}

void ConnectivityMatrix::CalculateMeanSignals( MatrixType& meanMatrix )
{
  // access input matrix
  PixelType* buffer = m_Input->GetPixelContainer()->GetBufferPointer();
  InputMatrixType inputMatrix( m_NumberOfSamples, m_NumberOfVoxels, buffer );

  // label vector
  LabelPixelType* labels = m_Labels->GetPixelContainer()->GetBufferPointer();

  // number of samples
  int sampleCount = m_NumberOfSamples - m_SkipInitial - m_SkipEnd;

  // mean matrix
  meanMatrix = MatrixType( sampleCount, m_NumberOfLabels );
  meanMatrix.fill( 0 );

  // keep track of number of voxels
  std::vector< int > labelSize;
  for( int i = 0; i < m_NumberOfLabels; ++i )
    {
    labelSize.push_back( 0 );
    }

  // add label voxels
  for( int voxel = 0; voxel < m_NumberOfVoxels; ++voxel )
    {
    if ( labels[ voxel ] < 1 || m_LabelSet.find( labels[ voxel ] ) == m_LabelSet.end() )
      {
      continue;
      }

    int column = m_LabelSet[ labels[ voxel ] ];
    labelSize[ column ]++;

    for( int row = m_SkipInitial; row < ( m_NumberOfSamples - m_SkipEnd ); ++row )
      {
      meanMatrix( row - m_SkipInitial, column ) += inputMatrix( row, voxel );
      }
    }

  m_NumberOfSamples = sampleCount;

  // average
  for( int i = 0; i < m_NumberOfLabels; ++i )
    {
    PixelType factor = 1.0 / static_cast< PixelType >( labelSize[ i ] );
    for( int j = 0; j < m_NumberOfSamples; ++j )
      {
      meanMatrix( j, i ) *= factor;
      }
    }
}

void ConnectivityMatrix::CalculatePCASignals( MatrixType& meanMatrix )
{
  // access input matrix
  PixelType* buffer = m_Input->GetPixelContainer()->GetBufferPointer();
  InputMatrixType inputMatrix( m_NumberOfSamples, m_NumberOfVoxels, buffer );

  // label vector
  LabelPixelType* labels = m_Labels->GetPixelContainer()->GetBufferPointer();

  // number of samples
  int sampleCount = m_NumberOfSamples - m_SkipInitial - m_SkipEnd;

  // mean matrix
  meanMatrix = MatrixType( sampleCount, m_NumberOfLabels );
  meanMatrix.fill( 0 );

  // keep track of number of voxels
  std::vector< int > labelSize;
  for( int i = 0; i < m_NumberOfLabels; ++i )
    {
    labelSize.push_back( 0 );
    }

  for( LabelSetType::iterator i = m_LabelSet.begin(); i != m_LabelSet.end(); ++i )
    {
    LabelPixelType label = i->first;
    int column = i->second;

    int numberOfSamples = 0;
    for( int voxel = 0; voxel < m_NumberOfVoxels; ++voxel )
      {
      if ( labels[ voxel ] != label )
        {
        continue;
        }

      ++numberOfSamples;
      }

    MatrixType sampleMatrix( sampleCount, numberOfSamples );

    int sample = 0;
    for( int voxel = 0; voxel < m_NumberOfVoxels; ++voxel )
      {
      if ( labels[ voxel ] != label )
        {
        continue;
        }

      for( int row = m_SkipInitial; row < ( m_NumberOfSamples - m_SkipEnd ); ++row )
        {
        sampleMatrix( row - m_SkipInitial, sample ) = inputMatrix( row, voxel );
        }

      ++sample;
      }

    PCA pca( sampleMatrix );
    pca.TransformPCA();
    VectorType signal = pca.GetComponent( 0 );
//
//    PixelType meanCC = 0;
//    PixelType maxCC = -1;
//    PixelType minCC = 1;
//
//    VectorType meanVector( sampleMatrix.rows() );
//    meanVector.fill( 0 );
//    for( unsigned int j = 0; j < sampleMatrix.cols(); ++j )
//      {
//      PixelType cc = rs::Common< PixelType >::CorrelationCoefficient( signal, sampleMatrix.get_column( j ) );
//      meanCC += cc;
//      maxCC = maxCC > cc ? maxCC : cc;
//      minCC = minCC < cc ? minCC : cc;
//      meanVector += sampleMatrix.get_column( j );
//      }
//
//    meanCC /= static_cast< PixelType >( sampleMatrix.cols() );
//    meanVector /= static_cast< PixelType >( sampleMatrix.cols() );
//
//    PixelType meanVectorCC = rs::Common< PixelType >::CorrelationCoefficient( signal, meanVector );
//    if ( meanVectorCC < 0 )
//      {
//      signal = signal * -1.;
//      }
//
    meanMatrix.set_column( column, signal );

//    std::cout << "Label " << label << ", column " << column << ", " << sampleMatrix.rows() << "x" << sampleMatrix.cols() << "; PCA CC=[" << minCC << ", " << meanCC << ", " << maxCC << "] / " << meanVectorCC << std::endl;
    }

  m_NumberOfSamples = sampleCount;
}

void ConnectivityMatrix::Run( bool doPCA )
{
  // initialize
  this->Initialize();

  // calculate for each label the ROI mean signal
  MatrixType meanMatrix;

  if ( doPCA )
    {
    this->CalculatePCASignals( meanMatrix );
    }
  else
    {
    this->CalculateMeanSignals( meanMatrix );
    }

  // output matrix
  PixelType* buffer = m_Output->GetPixelContainer()->GetBufferPointer();
  InputMatrixType outputMatrix( m_NumberOfLabels, m_NumberOfLabels, buffer );

  // calculate correlation coefficients between each ROI mean
  for( int i = 0; i < m_NumberOfLabels; ++i )
    {
    VectorType a = meanMatrix.get_column( i );
    for( int j = i + 1; j < m_NumberOfLabels; ++j )
      {
      VectorType b = meanMatrix.get_column( j );

      PixelType r = Common< PixelType >::CorrelationCoefficient( a, b );

      outputMatrix( i, j ) = r;
      outputMatrix( j, i ) = r;
      }
    }
}

ConnectivityMatrix::OutputType::Pointer ConnectivityMatrix::GetOutput()
{
  return m_Output;
}

ConnectivityMatrix::LabelSetType ConnectivityMatrix::GetLabelSet()
{
  return m_LabelSet;
}

} // end namespace rs

