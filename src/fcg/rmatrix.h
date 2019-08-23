#ifndef __fcgRMatrix_h__
#define __fcgRMatrix_h__
#include "itkImage.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "itkImageMaskSpatialObject.h"
#include "vnl/vnl_matrix_ref.h"
#include <vector>

namespace fcg
{

class RMatrix
{
public:
  typedef float PixelType;
  typedef float OutputPixelType;
  typedef unsigned char MaskPixelType;

  typedef itk::Image< PixelType, 4 > ImageType;
  typedef itk::Image< MaskPixelType, 3 > MaskImageType;
  typedef itk::ImageMaskSpatialObject< 3 > MaskType;

  typedef vnl_vector< PixelType > VectorType;
  typedef vnl_matrix_ref< PixelType > DataType;
  typedef vnl_matrix< PixelType > MatrixType;
  typedef itk::Image< PixelType, 3 > MapType;
  typedef vnl_matrix_ref< MaskPixelType > DataMaskType;
  typedef itk::Image< OutputPixelType, 2 > OutputImageType;
  typedef vnl_matrix_ref< OutputPixelType > OutputDataType;

  void Run( const std::string& inputFileName, const std::string& maskFileName, const std::string& outputFileName, bool transformToZ, int skipBegin = 50, int skipEnd = 50, float constant = 0 );

protected:
  ImageType::Pointer LoadImage( const std::string& filename );
  MaskType::Pointer LoadMask( const std::string& filename );
};

} // end namespace fcg

#endif /*__fcgRMatrix_h__*/
