#ifndef __rsConnectivityMatrix_h__
#define __rsConnectivityMatrix_h__
#include "itkImage.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include <map>

namespace rs
{

class ConnectivityMatrix
{
public:
  typedef float PixelType;
  typedef unsigned short LabelPixelType;
  typedef itk::Image< PixelType, 4 > InputType;
  typedef itk::Image< PixelType, 2 > OutputType;
  typedef itk::Image< LabelPixelType, 3 > LabelType;
  typedef vnl_matrix_ref< PixelType > InputMatrixType;
  typedef vnl_matrix< PixelType > MatrixType;
  typedef vnl_vector< PixelType > VectorType;
  typedef std::map< LabelPixelType, int > LabelSetType;

  ConnectivityMatrix();
  virtual ~ConnectivityMatrix();

  void SetInput( InputType::Pointer input );
  void SetLabels( LabelType::Pointer label );
  void SetLabelRange( LabelPixelType min, LabelPixelType max );

  void SetSkipTimePoints( int initial, int end );

  void Run( bool doPCA = false );

  OutputType::Pointer GetOutput();
  LabelSetType GetLabelSet();

protected:
  void Initialize();
  void CalculateMeanSignals( MatrixType& meanMatrix );
  void CalculatePCASignals( MatrixType& meanMatrix );

protected:
  int m_NumberOfLabels;
  int m_NumberOfVoxels;
  int m_NumberOfSamples;
  int m_SkipInitial;
  int m_SkipEnd;
  LabelPixelType m_LabelMin;
  LabelPixelType m_LabelMax;
  InputType::Pointer m_Input;
  LabelType::Pointer m_Labels;
  OutputType::Pointer m_Output;
  LabelSetType m_LabelSet;
};

} // end namespace rs

#endif /*__rsConnectivityMatrix_h__*/

