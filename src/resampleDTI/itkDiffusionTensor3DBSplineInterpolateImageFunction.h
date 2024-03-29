/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DBSplineInterpolateImageFunction.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DBSplineInterpolateImageFunction_h
#define __itkDiffusionTensor3DBSplineInterpolateImageFunction_h

#include "itkDiffusionTensor3DInterpolateImageFunctionReimplementation.h"
#include <itkBSplineInterpolateImageFunction.h>
#include <itkImage.h>

namespace itk
{
/**
 * \class DiffusionTensor3DSplineInterpolateImageFunction
 * 
 * Implementation of blockwise spline interpolation for diffusion tensor images
 */
template< class TData >
class DiffusionTensor3DBSplineInterpolateImageFunction :
  public DiffusionTensor3DInterpolateImageFunctionReimplementation< TData >
{
public:
  typedef TData DataType ;
  typedef DiffusionTensor3DBSplineInterpolateImageFunction Self ;
  typedef DiffusionTensor3DInterpolateImageFunctionReimplementation< DataType > Superclass ;
  typedef typename Superclass::ImageType ImageType ;
  typedef SmartPointer< Self > Pointer ;
  typedef SmartPointer< const Self > ConstPointer ;
  typedef BSplineInterpolateImageFunction< ImageType , double , double >
                                                 BSplineInterpolateFunction ;
  
  itkNewMacro( Self ) ;
  ///Get the Spline Order, supports 0th - 5th order splines. The default is a 1st order spline. 
  itkGetMacro( SplineOrder , int ) ;
  ///Set the Spline Order, supports 0th - 5th order splines. The default is a 1st order spline.
  itkSetMacro( SplineOrder , unsigned int ) ;

  protected:
  void AllocateInterpolator() ;
  DiffusionTensor3DBSplineInterpolateImageFunction() ;
  unsigned int m_SplineOrder ;
  typename BSplineInterpolateFunction::Pointer bSplineInterpolateFunction[ 6 ] ;
};

}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DBSplineInterpolateImageFunction.txx"
#endif

#endif
