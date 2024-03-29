/*=========================================================================

  Program:   Diffusion Applications
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/trunk/Applications/CLI/DiffusionApplications/ResampleDTI/itkMatrixExtended.h $
  Language:  C++
  Date:      $Date: 2008-11-25 20:23:08 +0100 (Tue, 25 Nov 2008) $
  Version:   $Revision: 7976 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/
#ifndef __itkMatrixExtended_h
#define __itkMatrixExtended_h

#include "itkMatrix.h"

namespace itk
{
/** \class MatrixExtended
 *
 * Matrices of this class can be casted blockwise 
 */    
template< class T , unsigned int NRows = 3 , unsigned int NColumns = 3 >
class MatrixExtended : public Matrix< T , NRows , NColumns >
{
public:
  typedef MatrixExtended Self ;
  typedef Matrix< T , NRows , NColumns > Superclass ;
  typedef T ValueType ;
  typedef T ComponentType ;
  typedef vnl_matrix_fixed< T , NRows , NColumns > InternalMatrixType ;
  MatrixExtended() ;
  MatrixExtended( const InternalMatrixType &matrix ) ;
  MatrixExtended( const Self &matrix ) ;
  MatrixExtended( const Superclass &matrix ) ;
  ///Cast the matrix
  template< class C , unsigned int NRowsC , unsigned int NColumnsC >
  operator MatrixExtended < C , NRowsC , NColumnsC > const() ;
  Self operator=( const Self & matrix ) ;
  Self operator=( const Superclass & matrix ) ;
  Self operator=( const InternalMatrixType & matrix ) ;
};

}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMatrixExtended.txx"
#endif

#endif
