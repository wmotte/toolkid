#include "rsPCA.h"

namespace rs
{

PCA::PCA( const MatrixType& x ) : m_X( x )
{
}

void PCA::TransformPCA()
{
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

    m_Variance.clear();
    m_CumulativeVariance.clear();

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
  m_Components = m_L.size();

  for( unsigned int i = 0; i < m_Components; ++i )
    {
    PixelType variance = m_L( i ) / sum;
    cumulative += variance;
    m_CumulativeVariance.push_back( cumulative );
    m_Variance.push_back( variance );
    }

  // extract components
  m_D = svd.U();
}

PCA::VectorType PCA::GetComponent( unsigned int i ) const
{
  return m_D.get_column( i );
}

unsigned int PCA::GetNumberOfComponents() const
{
  return m_D.cols();
}


PCA::PixelType PCA::GetVariance( unsigned int i ) const
{
  return m_Variance[ i ];
}

PCA::PixelType PCA::GetCumulativeVariance( unsigned int i ) const
{
  return m_CumulativeVariance[ i ];
}

} // end namespace rs
