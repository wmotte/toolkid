/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DNonRigidTransform.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DNonRigidTransform_h
#define __itkDiffusionTensor3DNonRigidTransform_h

#include "itkDiffusionTensor3DTransform.h"
#include "itkDiffusionTensor3DFSAffineTransform.h"
#include "itkDiffusionTensor3DPPDAffineTransform.h"
#include <itkTransform.h>

namespace itk
{


template< class TData >
class DiffusionTensor3DNonRigidTransform :
  public DiffusionTensor3DTransform< TData >
{
public:
  typedef TData DataType ; 
  typedef DiffusionTensor3DNonRigidTransform Self ;
  typedef DiffusionTensor3DTransform< DataType > Superclass ;
  typedef typename Superclass::TensorDataType TensorDataType ;
  typedef typename Superclass::MatrixTransformType MatrixTransformType ;
  typedef typename Superclass::PointType PointType ;
  typedef SmartPointer< Self > Pointer ;
  typedef SmartPointer< const Self > ConstPointer ;
  typedef Transform< double, 3, 3 > TransformType ;
  typedef itk::DiffusionTensor3DPPDAffineTransform< DataType > PPDAffineTransformType ;
  typedef itk::DiffusionTensor3DFSAffineTransform< DataType > FSAffineTransformType ;
  typedef itk::DiffusionTensor3DAffineTransform< DataType > AffineTransform;
  //SmartPointer
  itkNewMacro( Self ) ;
  ///Set the transform
  itkSetObjectMacro( Transform , TransformType ) ;
  ///Evaluate the position of the transformed tensor in the output image
  PointType EvaluateTensorPosition( const PointType &point ) ;
  ///Evaluate the transformed tensor
  TensorDataType EvaluateTransformedTensor( TensorDataType &tensor , PointType &outputPosition ) ;
  void SetAffineTransformType(typename AffineTransform::Pointer transform);

protected:
  unsigned long latestTime ;
  typename TransformType::Pointer m_Transform ;
  typename AffineTransform::Pointer m_Affine;
};


}//end of itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DNonRigidTransform.txx"
#endif

#endif

