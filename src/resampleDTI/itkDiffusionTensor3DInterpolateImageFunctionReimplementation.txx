/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkDiffusionTensor3DInterpolateImageFunctionReimplementation.txx $
  Language:  C++
  Date:      $Date: 2009-07-24 17:09:59 +0200 (Fri, 24 Jul 2009) $
  Version:   $Revision: 10016 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __DiffusionTensor3DInterpolateImageFunctionReimplementation_txx
#define __DiffusionTensor3DInterpolateImageFunctionReimplementation_txx

#include "itkDiffusionTensor3DInterpolateImageFunctionReimplementation.h"


namespace itk
{
    

template< class TData >
void
DiffusionTensor3DInterpolateImageFunctionReimplementation< TData >
::PreCompute()
{
  typename DiffusionImageType::RegionType region ;
  region = this->m_InputImage->GetLargestPossibleRegion() ;
  typename DiffusionImageType::SizeType size = region.GetSize() ;
  typename DiffusionImageType::PointType origin = this->m_InputImage->GetOrigin() ;
  typename DiffusionImageType::DirectionType direction = this->m_InputImage->GetDirection() ;
  typename DiffusionImageType::SpacingType spacing = this->m_InputImage->GetSpacing() ; 
  IteratorDiffusionImageType it( this->m_InputImage , size ) ;
  std::vector< IteratorImageType > out ;
  for( int i = 0 ; i < 6 ; i++ ) 
    {
    m_Image[ i ] = ImageType::New() ;
    m_Image[ i ]->SetRegions( size ) ;
    m_Image[ i ]->SetOrigin( origin ) ;
    m_Image[ i ]->SetDirection( direction ) ;
    m_Image[ i ]->SetSpacing( spacing ) ;
    m_Image[ i ]->Allocate() ;
    IteratorImageType outtemp( m_Image[ i ] , m_Image[ i ]->GetLargestPossibleRegion() ) ;
    out.push_back( outtemp ) ;
    }
  for( int i = 0 ; i < 6 ; i++ )
    { out[ i ].GoToBegin() ; }
  for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
    {
    TensorDataType tensor = it.Get() ;
    for( int i = 0 ; i < 6 ; i++ )
      {
      out[ i ].Set( tensor[ i ] ) ;
      ++out[ i ] ;
      }
    }
  AllocateInterpolator() ;
  for( int i = 0 ; i < 6 ; i++ )
    {
    interpol[ i ]->SetInputImage( m_Image[ i ] ) ;
    }
  //Compute position of the lower and superior corner of the image    
  this->PreComputeCorners() ; 
  this->latestTime = Object::GetMTime() ;
}
    
template< class TData >
typename DiffusionTensor3DInterpolateImageFunctionReimplementation< TData >
::TensorDataType 
DiffusionTensor3DInterpolateImageFunctionReimplementation< TData >
::Evaluate( const PointType &point )
{
  TensorDataType pixelValue( NumericTraits< DataType >::Zero ) ;
  if( this->m_InputImage.IsNotNull() )
    {
    if( this->latestTime< Object::GetMTime() )
      { 
      this->P->Down() ;
      if( this->latestTime< Object::GetMTime() )
        {
        PreCompute() ;
        }
      this->P->Up() ;
      }
    bool ok=1 ;
    if( !interpol[0]->IsInsideBuffer( point) )
      {
      ok = false;
      }
    /***  Replaced by the code above to account for directional vectors 
    for( int i = 0 ; i < 3 ; i++ )
      {
      if( point[ i ] >this->m_End[ i ] || point[ i ]< this->m_Origin[ i ] )
        { ok = 0 ; }
      }
    **/
    if( ok )
      {
      for( int i = 0 ; i < 6 ; i++ )
        {
        pixelValue[ i ] = ( DataType ) interpol[ i ]->Evaluate( point ) ;
        }
      }
    }
  else
    {
    itkExceptionMacro( << "No InputImage Set" ) ;
    }  
  return pixelValue ;
}

}//end itk namespace

#endif
