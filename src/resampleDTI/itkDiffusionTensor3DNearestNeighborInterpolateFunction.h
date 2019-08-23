/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DNearestNeighborInterpolateFunction.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DNearestNeighborInterpolateFunction_h
#define __itkDiffusionTensor3DNearestNeighborInterpolateFunction_h


#include "itkDiffusionTensor3DInterpolateImageFunction.h"

namespace itk
{

/** \class DiffusionTensor3DInterpolateImageFunction
 * 
 * Implementation of the nearest neighborhood interpolation for diffusion tensor images
 */

template< class TData >
class DiffusionTensor3DNearestNeighborInterpolateFunction :
  public DiffusionTensor3DInterpolateImageFunction< TData >
{
public:
  typedef TData DataType ;
  typedef DiffusionTensor3DNearestNeighborInterpolateFunction Self ;
  typedef DiffusionTensor3DInterpolateImageFunction< DataType > Superclass ;
  typedef typename Superclass::TensorDataType TensorDataType ;
  typedef typename Superclass::DiffusionImageType DiffusionImageType ;
  typedef typename Superclass::PointType PointType ;
  typedef SmartPointer< Self > Pointer ;
  typedef SmartPointer< const Self > ConstPointer ;

  itkNewMacro( Self ) ;
  ///Evaluate the value of a tensor at a given position
  TensorDataType Evaluate( const PointType &point ) ;

};

}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DNearestNeighborInterpolateFunction.txx"
#endif

#endif
