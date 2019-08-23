#include "tspca.h"

PCA::PCA( unsigned int components, float variance ) : m_Components( components ), m_Variance( variance )
{
}

void PCA::ReadImage( const std::string& filename )
{
  std::cout << "Read " << filename << std::endl;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();

  m_Image = reader->GetOutput();
}

void PCA::ReadMask( const std::string& filename )
{
  std::cout << "Read " << filename << std::endl;
  typedef itk::ImageFileReader< MaskType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();

  m_Mask = reader->GetOutput();
}

void PCA::ImageToTimeSeries()
{
  MaskType::Pointer mask = MaskType::New();
  MaskType::RegionType maskRegion;
  MaskType::SizeType maskSize;
  MaskType::IndexType maskIndex;
  MaskType::SpacingType maskSpacing;
  MaskType::DirectionType maskDirection;
  MaskType::PointType maskOrigin;
  ImageType::RegionType imageRegion = m_Image->GetLargestPossibleRegion();
  ImageType::SizeType imageSize = imageRegion.GetSize();
  ImageType::DirectionType imageDirection = m_Image->GetDirection();
  ImageType::PointType imageOrigin = m_Image->GetOrigin();
  ImageType::SpacingType imageSpacing = m_Image->GetSpacing();
  for( int i = 0; i < 3; ++i )
    {
    maskSize[ i ] = imageSize[ i ];
    maskSpacing[ i ] = imageSpacing[ i ];
    maskOrigin[ i ] = imageOrigin[ i ];
    maskIndex[ i ] = 0;

    for( int j = 0; j < 3; ++j )
      {
      maskDirection( i, j ) = imageDirection( i, j );
      }
    }
  maskRegion.SetSize( maskSize );
  maskRegion.SetIndex( maskIndex );
  mask->SetSpacing( maskSpacing );
  mask->SetOrigin( maskOrigin );
  mask->SetDirection( maskDirection );
  mask->SetRegions( maskRegion );
  mask->Allocate();
  mask->FillBuffer( 0 );
  MaskType::RegionType sourceRegion = m_Mask->GetLargestPossibleRegion();

  typedef itk::LinearInterpolateImageFunction< MaskType, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( m_Mask );

  int samples = 0;
  itk::ImageRegionIteratorWithIndex< MaskType > itMask( mask, maskRegion );
  for( ; !itMask.IsAtEnd(); ++itMask )
    {
    MaskType::PointType maskPoint;
    maskIndex = itMask.GetIndex();
    mask->TransformIndexToPhysicalPoint( maskIndex, maskPoint );
    m_Mask->TransformPhysicalPointToIndex( maskPoint, maskIndex );

//      if ( sourceRegion.IsInside( maskIndex ) && m_Mask->GetPixel( maskIndex ) > 0 )
    if ( interpolator->IsInsideBuffer( maskPoint ) && interpolator->Evaluate( maskPoint ) > 0 )
      {
      itMask.Set( 255 );
      ++samples;
      }
    }

  m_Mask = mask;

  ImageType::RegionType region = m_Image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();

  int timepoints = size[ 3 ];

  std::cout << timepoints << " x " << samples << " timepoints x samples" << std::endl;

  m_X = MatrixType( timepoints, samples );

  ImageType::IndexType index = region.GetIndex();
  for( int i = 0; i < timepoints; ++i )
    {
    index[ 3 ] = i;
    size[ 3 ] = 0;
    region.SetSize( size );
    region.SetIndex( index );
    itk::ImageRegionConstIterator< ImageType > itImage( m_Image, region );

    int sample = 0;
    for( itMask.GoToBegin(); !itMask.IsAtEnd(); ++itMask, ++itImage )
      {
      if ( itMask.Value() > 0 )
        {
        m_X( i, sample++ ) = itImage.Value();
        }
      }
    }
}

void PCA::TransformPCA()
{
  std::cout << "Subtract mean" << std::endl;

  // not really?
//  m_Mean = VectorType( m_X.cols() );
//  for( unsigned int j = 0; j < m_X.cols(); ++j )
//    {
//    VectorType column = m_X.get_column( j );
//    PixelType mean = column.mean();
//    for( unsigned int i = 0; i < m_X.rows(); ++i )
//      {
//      m_X( i, j ) -= mean;
//      }
//    m_Mean( j ) = mean;
//    }

    m_Mean = VectorType( m_X.rows() );
    m_Mean.fill( 0 );
    for( unsigned int i = 0; i < m_X.cols(); ++i )
      {
      m_Mean += m_X.get_column( i );
      }

    m_Mean /= static_cast< PixelType >( m_X.cols() );

    for( unsigned int i = 0; i < m_X.rows(); ++i )
      {
      for( unsigned int j = 0; j < m_X.cols(); ++j )
        {
        m_X( i, j ) -= m_Mean( i );
        }
      }

//    std::cout << "Normalize time courses by standard deviation" << std::endl;
//    for( int i = 0; i < m_X.cols(); ++i )
//      {
//      VectorType column = m_X.get_column( i );
//      PixelType mean = column.mean();
//
//      PixelType sumsq = 0;
//      for( int j = 0; j < m_X.rows(); ++j )
//        {
//        PixelType entry = column( j ) - mean;
//        sumsq += entry * entry;
//        }
//
//      PixelType std = vcl_sqrt( sumsq / static_cast< PixelType >( m_X.rows() - 1 ) );
//
//      for( int j = 0; j < m_X.rows(); ++j )
//        {
//        m_X( j, i ) /= std;
//        }
//      }

  std::cout << "Calculate covariance matrix" << std::endl;
  MatrixType cov = MatrixType( m_X.rows(), m_X.rows() );
  cov.fill( 0 );

  // X * X', compute upper triangle
  for( unsigned int i = 0; i < m_X.rows(); ++i )
    {
    for( unsigned int j = i; j < m_X.rows(); ++j )
      {
      for( unsigned int k = 0; k < m_X.cols(); ++k )
        {
        cov( i, j ) += m_X( i, k ) * m_X( j, k );
        }
      }
    }

  // copy lower triangle
  for( unsigned int j = 0; j < m_X.rows(); ++j )
    {
    for( unsigned int i = j + 1; i < m_X.rows(); ++i )
      {
      cov( i, j ) = cov( j, i );
      }
    }

  // normalize
  cov /= static_cast< PixelType >( m_X.rows() );

  std::cout << "Singular value decomposition" << std::endl;
  vnl_svd< PixelType > svd( cov );

  // get diagonal elements
  m_L = svd.W().diagonal();

  // count variance
  PixelType sum = 0;
  for( unsigned int i = 0; i < m_L.size(); ++i )
    {
    sum += m_L( i );
    }

  // find cut-off variance
  PixelType cumulative = 0;
  int numberOfVarianceComponents = 0;
  for( unsigned int i = 0; i < m_L.size() && cumulative < m_Variance; ++i, ++numberOfVarianceComponents )
    {
    cumulative += m_L( i ) / sum;
    }

  // components
  unsigned int components = ( m_Components < 1 || m_Components > m_X.cols() ) ? m_X.cols() : m_Components;
  if ( m_Variance > 0 )
    {
    components = numberOfVarianceComponents;
    std::cout << "Extract " << ( cumulative * 100. ) << "%" << std::endl;
    }

  m_Components = components;
  std::cout << "Extract " << components << " components" << std::endl;

  // extract sub-matrix
  m_D = svd.U().extract( m_X.rows(), components );

  std::cout << "Solve for regression coefficients" << std::endl;
  // D*b = X -> obtain b
  vnl_svd< PixelType > svdD( m_D );
  m_b = svdD.solve( m_X );
}

void PCA::WriteReconstructedImage( const std::string& filename )
{
  std::cout << "Write reconstructed image to " << filename << std::endl;
  ImagePointer output = ImageType::New();
  output->CopyInformation( m_Image );
  itk::ImageRegionConstIterator< MaskType > itMask( m_Mask, m_Mask->GetLargestPossibleRegion() );
  ImageType::RegionType region = m_Image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType index = region.GetIndex();

  MatrixType X = m_D * m_b;

  output->SetRegions( region );
  output->Allocate();
  output->FillBuffer( 0 );

  int timepoints = size[ 3 ];

  for( int i = 0; i < timepoints; ++i )
    {
    size[ 3 ] = 0;
    index[ 3 ] = i;

    region.SetSize( size );
    region.SetIndex( index );

    itk::ImageRegionIterator< ImageType > it( output, region );
    int sample = 0;
    for( itMask.GoToBegin(); !itMask.IsAtEnd(); ++it, ++itMask )
      {
      if ( itMask.Value() > 0 )
        {
        it.Set( m_Mean( i ) + X( i, sample++ ) );
        }
      else
        {
        it.Set( 0 );
        }
      }
    }

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( output );
  writer->Update();
}



void PCA::WriteDesignMatrix( const std::string& filename )
{
  std::cout << "Write design matrix to " << filename << std::endl;
  std::ofstream out( filename.c_str() );
  for( unsigned int i = 0; i < m_D.rows(); ++i )
    {
    for( unsigned int j = 0; j < m_D.cols(); ++j )
      {
      out << ( j == 0 ? "" : "\t" ) << m_D( i, j );
      }

    out << std::endl;
    }
}

void PCA::WriteRegressionCoefficients( const std::string& filename )
{
  std::cout << "Write coefficient image to " << filename << std::endl;
  ImagePointer output = ImageType::New();
  output->CopyInformation( m_Image );
  itk::ImageRegionConstIterator< MaskType > itMask( m_Mask, m_Mask->GetLargestPossibleRegion() );
  ImageType::RegionType region = m_Image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType index = region.GetIndex();

//    int timepoints = size[ 3 ];
  int timepoints = m_b.rows();

  size[ 3 ] = timepoints;
  region.SetSize( size );

  output->SetRegions( region );
  output->Allocate();
  output->FillBuffer( 0 );


  for( int i = 0; i < timepoints; ++i )
    {
    size[ 3 ] = 0;
    index[ 3 ] = i;

    region.SetSize( size );
    region.SetIndex( index );

    itk::ImageRegionIterator< ImageType > it( output, region );
    int sample = 0;
    for( itMask.GoToBegin(); !itMask.IsAtEnd(); ++it, ++itMask )
      {
      if ( itMask.Value() > 0 )
        {
        it.Set( m_b( i, sample++ ) );
        }
      else
        {
        it.Set( 0 );
        }
      }
    }

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( output );
  writer->Update();
}

void PCA::WriteCoefficientMatrix( const std::string& filename )
{
  std::cout << "Write coefficient matrix to " << filename << std::endl;
  std::ofstream out( filename.c_str() );
  for( unsigned int i = 0; i < m_b.rows(); ++i )
    {
    for( unsigned int j = 0; j < m_b.cols(); ++j )
      {
      out << ( j == 0 ? "" : "\t" ) << m_b( i, j );
      }

    out << std::endl;
    }
}

void PCA::WriteCorrelationImage( const std::string& filename )
{
  std::cout << "Write correlation image to " << filename << std::endl;
  ImageType::RegionType region = m_Image->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType index = region.GetIndex();

  unsigned int components = m_D.cols();
  size[ 3 ] = components;
  region.SetSize( size );

  ImagePointer output = ImageType::New();
  output->CopyInformation( m_Image );
  output->SetRegions( region );
  output->Allocate();
  output->FillBuffer( 0 );

  for( unsigned int i = 0; i < components; ++i )
    {
    VectorType tc = m_D.get_column( i );

    size[ 3 ] = 0;
    index[ 3 ] = i;

    region.SetSize( size );
    region.SetIndex( index );

    itk::ImageRegionIterator< ImageType > itOut( output, region );
    itk::ImageRegionConstIterator< MaskType > itMask( m_Mask, m_Mask->GetLargestPossibleRegion() );

    int sample = 0;
    for( ; !itMask.IsAtEnd(); ++itMask, ++itOut )
      {
      if ( itMask.Value() > 0 )
        {
        VectorType pattern = m_X.get_column( sample++ );
        PixelType correlation = CorrelationCoefficient( tc, pattern );
        itOut.Set( correlation );
        }
      }
    }

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( output );
  writer->Update();
}

PCA::PixelType PCA::CorrelationCoefficient( const VectorType& a, const VectorType& b )
{
  // cov( a, b ) / ( std( a ) * std( b ) )

  // cov( a, b ):
  // ( a - mean )' * ( b - mean ) / ( n - 1 )

  VectorType c = a - a.mean();
  VectorType d = b - b.mean();

  PixelType cov = dot_product< PixelType >( c, d ) / static_cast< PixelType >( c.size() - 1 );
  PixelType sumC = 0;
  PixelType sumD = 0;
  for( unsigned int i = 0; i < c.size(); ++i )
    {
    sumC += c( i ) * c( i );
    sumD += d( i ) * d( i );
    }

  PixelType stdC = vcl_sqrt( sumC / static_cast< PixelType >( c.size() - 1 ) );
  PixelType stdD = vcl_sqrt( sumD / static_cast< PixelType >( d.size() - 1 ) );

  return ( cov / ( stdC * stdD ) );
}

void PCA::WriteEigenvalues( const std::string& filename )
{
  std::cout << "Write explained variances to " << filename << std::endl;

  PixelType sum = 0;
  for( unsigned int i = 0; i < m_L.size(); ++i )
    {
    sum += m_L( i );
    }

  std::ofstream out( filename.c_str() );
  PixelType cumulative = 0;
  for( unsigned int i = 0; i < m_L.size(); ++i )
    {
    PixelType partial = ( m_L( i ) / sum );
    cumulative += partial;
    out << m_L( i ) << "\t" << partial << "\t" << cumulative << "\t" << ( 1 - cumulative ) << std::endl;
    }
}
