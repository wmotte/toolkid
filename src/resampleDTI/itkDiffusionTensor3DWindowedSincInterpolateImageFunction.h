/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DWindowedSincInterpolateImageFunction.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DWindowedSincInterpolateImageFunction_h
#define __itkDiffusionTensor3DWindowedSincInterpolateImageFunction_h

#include "itkDiffusionTensor3DInterpolateImageFunctionReimplementation.h"
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkImage.h>

namespace itk
{
/**
 * \class DiffusionTensor3DWindowedSincInterpolateImageFunction
 * 
 * Implementation of windowed sinc blockwise interpolation for diffusion tensor images
 */
template< class TData,
         unsigned int VRadius ,
         class TWindowFunction = Function::HammingWindowFunction< VRadius > ,
         class TBoundaryCondition = ConstantBoundaryCondition< OrientedImage< TData , 3 > > >
class DiffusionTensor3DWindowedSincInterpolateImageFunction
 : public DiffusionTensor3DInterpolateImageFunctionReimplementation< TData >
{
public:
  typedef TData DataType ;
  typedef DiffusionTensor3DWindowedSincInterpolateImageFunction Self ;
  typedef DiffusionTensor3DInterpolateImageFunctionReimplementation< DataType > Superclass ;
  typedef typename Superclass::ImageType ImageType ;
  typedef SmartPointer< Self > Pointer ;
  typedef SmartPointer< const Self > ConstPointer ;
  typedef WindowedSincInterpolateImageFunction< ImageType ,
                                    VRadius , TWindowFunction,
                                    TBoundaryCondition , double > WindowedSincInterpolateImageFunctionType ;
  
  itkNewMacro(Self);

protected:
  void AllocateInterpolator() ;
  typename WindowedSincInterpolateImageFunctionType::Pointer windowedSincInterpolator[ 6 ] ;
};

}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DWindowedSincInterpolateImageFunction.txx"
#endif

#endif
