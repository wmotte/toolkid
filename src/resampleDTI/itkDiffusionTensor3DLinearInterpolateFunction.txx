/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DLinearInterpolateFunction.txx $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkDiffusionTensor3DLinearInterpolateFunction_txx
#define __itkDiffusionTensor3DLinearInterpolateFunction_txx

#include "itkDiffusionTensor3DLinearInterpolateFunction.h"

namespace itk
{

    
template< class TData >
void
DiffusionTensor3DLinearInterpolateFunction< TData >
::AllocateInterpolator()
{
  for( int i = 0 ; i < 6 ; i++ )
    {
    linearInterpolator[ i ] = LinearInterpolateImageFunctionType::New() ;
    this->interpol[ i ] = linearInterpolator[ i ] ;
    }
}

}//end itk namespace

#endif
