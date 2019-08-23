#ifndef __tspca_h__
#define __tspca_h__
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_svd.h"
#include "vcl_cmath.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "tkdCmdParser.h"

class PCA
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 4 > ImageType;
  typedef itk::Image< unsigned char, 3 > MaskType;
  typedef ImageType::Pointer ImagePointer;
  typedef MaskType::Pointer MaskPointer;
  typedef vnl_matrix< float > MatrixType;
  typedef vnl_vector< float > VectorType;

  PCA( unsigned int components = 0, float variance = 0 );
  void ReadImage( const std::string& filename );
  void ReadMask( const std::string& filename );
  void ImageToTimeSeries();
  void TransformPCA();
  void WriteReconstructedImage( const std::string& filename );
  void WriteDesignMatrix( const std::string& filename );
  void WriteRegressionCoefficients( const std::string& filename );
  void WriteCoefficientMatrix( const std::string& filename );
  void WriteCorrelationImage( const std::string& filename );
  PixelType CorrelationCoefficient( const VectorType& a, const VectorType& b );
  void WriteEigenvalues( const std::string& filename );

protected:
  ImagePointer m_Image;
  MaskPointer m_Mask;
  MatrixType m_X;
  MatrixType m_D;
  VectorType m_Mean;
  MatrixType m_b;
  VectorType m_L;
  unsigned int m_Components;
  float m_Variance;
};

#endif /*__tspca_h__*/
