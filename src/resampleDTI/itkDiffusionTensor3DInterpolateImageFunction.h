/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DInterpolateImageFunction.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DInterpolateImageFunction_h
#define __itkDiffusionTensor3DInterpolateImageFunction_h

#include <itkObject.h>
#include "itkDiffusionTensor3D.h"
#include <itkOrientedImage.h>
#include <itkPoint.h>
#include <itkSemaphore.h>
#include <itkNumericTraits.h>

namespace itk
{
/**
 * \class DiffusionTensor3DInterpolateImageFunction
 * 
 * Virtual class to implement diffusion tensor interpolation classes 
 * 
 */
template< class TData >
class DiffusionTensor3DInterpolateImageFunction : public Object
{
public :
  typedef TData TensorType ;
  typedef DiffusionTensor3DInterpolateImageFunction Self ;
  typedef DiffusionTensor3D< TensorType > TensorDataType ;
  typedef OrientedImage< TensorDataType , 3 > DiffusionImageType ;
  typedef typename DiffusionImageType::Pointer DiffusionImageTypePointer ;
  typedef Point< double , 3 > PointType ;
  typedef SmartPointer< Self > Pointer ;
  typedef SmartPointer< const Self > ConstPointer ;

  ///Set the input image
  itkSetObjectMacro( InputImage , DiffusionImageType ) ;
  ///Evaluate the tensor value at a given position
  virtual TensorDataType Evaluate( const PointType &point ) = 0 ;

protected:
  DiffusionTensor3DInterpolateImageFunction() ;
  DiffusionImageTypePointer m_InputImage ;
  Semaphore::Pointer P ;
  PointType m_Origin ;
  PointType m_End ;
  unsigned long latestTime ;
  void PreComputeCorners() ;
};

}//end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DInterpolateImageFunction.txx"
#endif

#endif
