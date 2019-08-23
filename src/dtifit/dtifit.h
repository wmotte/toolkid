#ifndef __dtifit_h___
#define __dtifit_h___
#include "itkImage.h"
#include "itkVectorImage.h"
#include "vnl/vnl_matrix.h"
#include "itkAutoPointer.h"
#include "itkRGBPixel.h"
#include "itkDiffusionTensor3D.h"

class Procparser;

namespace dtifit
{

class DTIFit : public itk::LightObject
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 4 > InputImageType;
  typedef vnl_matrix< double > GradientTableType;
  typedef itk::Image< PixelType, 3 > OutputImageType;
  typedef itk::Image< PixelType, 4 > OutputVectorImageType;
  typedef itk::DiffusionTensor3D< double > TensorPixelType;
  typedef itk::Image< TensorPixelType, 3 > TensorImageType;
  typedef itk::RGBPixel< unsigned char > RGBPixelType;
  typedef itk::Image< RGBPixelType, 3 > RGBImageType;

  typedef DTIFit Self;
  typedef itk::LightObject Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( DTIFit, LightObject );

  void SetBValue( double bvalue );
  void SetGradientTable( const GradientTableType& table );
  void SetImage( InputImageType::Pointer input );
  void Fit();

  OutputImageType::Pointer GetFA();
  OutputImageType::Pointer GetRA();
  OutputImageType::Pointer GetTrace();

  OutputImageType::Pointer GetL1();
  OutputImageType::Pointer GetL2();
  OutputImageType::Pointer GetL3();

  OutputVectorImageType::Pointer GetV1();
  OutputVectorImageType::Pointer GetV2();
  OutputVectorImageType::Pointer GetV3();

  TensorImageType::Pointer GetTensor();
  RGBImageType::Pointer GetColorFA();

  void Run( InputImageType::Pointer input, const Procparser& pp, const std::string& outputBase, const std::string& outputExtension = "nii.gz",  bool mirrorGradientTable = false );

protected:
  DTIFit();
  ~DTIFit();

private:
  struct Impl;
  itk::AutoPointer< Impl > m_Impl;

  DTIFit( const Self& );
  void operator=( const Self& );
};

} // end namespace dtifit

#endif /*__dtifit_h___*/
