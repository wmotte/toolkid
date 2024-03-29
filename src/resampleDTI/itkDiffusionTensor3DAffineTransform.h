/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DAffineTransform.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DAffineTransform_h
#define __itkDiffusionTensor3DAffineTransform_h


#include "itkDiffusionTensor3DMatrix3x3Transform.h"
#include <itkAffineTransform.h>


namespace itk
{
/**
 * \class DiffusionTensor3DAffineTransform
 * 
 * 
 * This is an abstract class to implement affine transformation for diffusion tensors.
 */
template< class TData >
class DiffusionTensor3DAffineTransform :
  public DiffusionTensor3DMatrix3x3Transform< TData >
{
public:
  typedef TData DataType ;
  typedef DiffusionTensor3DAffineTransform Self ;
  typedef DiffusionTensor3DMatrix3x3Transform< DataType > Superclass ;
  typedef typename Superclass::PointType PointType ;
  typedef typename Superclass::TensorDataType TensorDataType ;
  typedef typename Superclass::MatrixDataType MatrixDataType ;
  typedef typename Superclass::MatrixTransformType MatrixTransformType ;
  typedef typename Superclass::InternalTensorDataType InternalTensorDataType ;
  typedef typename Superclass::InternalMatrixDataType InternalMatrixDataType ;
  typedef typename Superclass::InternalMatrixTransformType InternalMatrixTransformType ;
  typedef SmartPointer< Self > Pointer ;
  typedef SmartPointer< const Self > ConstPointer ;
  typedef MatrixExtended< double , 4 , 4 > MatrixTransform4x4Type ;
  typedef AffineTransform< double , 3 > AffineTransformType ;
  typedef typename Superclass::VectorType VectorType;
  /** Set the transformation matrix from an itk::AffineTransform< double , 3 > object
  */
  void SetTransform( typename AffineTransformType::Pointer transform ) ;
  /** Set the transformation matrix directly from a 4x4 matrix
  * The 4th row is ignored and considered to be 0,0,0,1
  */
  void SetMatrix4x4( MatrixTransform4x4Type matrix ) ;

};

}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DAffineTransform.txx"
#endif

#endif
