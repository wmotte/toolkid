/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DBSplineInterpolateImageFunction.txx $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DBSplineInterpolateImageFunction_txx
#define __itkDiffusionTensor3DBSplineInterpolateImageFunction_txx

#include "itkDiffusionTensor3DBSplineInterpolateImageFunction.h"

namespace itk
{
    
template< class TData >
DiffusionTensor3DBSplineInterpolateImageFunction< TData >
::DiffusionTensor3DBSplineInterpolateImageFunction()
{
  m_SplineOrder = 1 ;
}    
    
template< class TData >
void
DiffusionTensor3DBSplineInterpolateImageFunction< TData >
::AllocateInterpolator()
{
  for( int i = 0 ; i < 6 ; i++ )
    {
    bSplineInterpolateFunction[ i ] = BSplineInterpolateFunction::New() ;
    bSplineInterpolateFunction[ i ]->SetSplineOrder( m_SplineOrder ) ;
    this->interpol[ i ] = bSplineInterpolateFunction[ i ] ;
    }
}


}//end itk namespace

#endif

