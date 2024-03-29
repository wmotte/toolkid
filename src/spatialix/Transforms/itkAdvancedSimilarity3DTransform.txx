/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdvancedSimilarity3DTransform.txx,v $
  Language:  C++
  Date:      $Date: 2007-11-27 16:04:48 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAdvancedSimilarity3DTransform_txx
#define _itkAdvancedSimilarity3DTransform_txx

#include "itkAdvancedSimilarity3DTransform.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_det.h"


namespace itk
{

// Constructor with default arguments
template <class TScalarType>
AdvancedSimilarity3DTransform<TScalarType>
::AdvancedSimilarity3DTransform() :
  Superclass(OutputSpaceDimension, ParametersDimension)
{
  m_Scale = 1.0;
  this->PrecomputeJacobianOfSpatialJacobian();
}


// Constructor with arguments
template<class TScalarType>
AdvancedSimilarity3DTransform<TScalarType>::
AdvancedSimilarity3DTransform( unsigned int outputSpaceDim, 
                        unsigned int paramDim) :
  Superclass(outputSpaceDim,paramDim)
{
}
 

// Constructor with arguments
template<class TScalarType>
AdvancedSimilarity3DTransform<TScalarType>::
AdvancedSimilarity3DTransform( const MatrixType & matrix,
                        const OutputVectorType & offset) :
  Superclass(matrix,offset)
{
}
 

// Set the scale factor
template <class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>
::SetScale( ScaleType scale )
{
   m_Scale = scale;
   this->ComputeMatrix();
}


// Directly set the matrix
template<class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>
::SetMatrix( const MatrixType & matrix )
{
  //
  // Since the matrix should be an orthogonal matrix
  // multiplied by the scale factor, then its determinant
  // must be equal to the cube of the scale factor.
  //
  double det = vnl_det( matrix.GetVnlMatrix() );

  if( det == 0.0 )
    {
    itkExceptionMacro( << "Attempting to set a matrix with a zero determinant" );
    }

  //
  // A negative scale is not acceptable
  // It will imply a reflection of the coordinate system.
  //

  double s = vnl_math_cuberoot( det );

  //
  // A negative scale is not acceptable
  // It will imply a reflection of the coordinate system.
  //
  if( s <= 0.0 )
    {
    itkExceptionMacro( << "Attempting to set a matrix with a negative trace" );
    }

  MatrixType testForOrthogonal = matrix;
  testForOrthogonal /= s;

  const double tolerance = 1e-10;
  if( !this->MatrixIsOrthogonal( testForOrthogonal, tolerance ) ) 
    {    
    itkExceptionMacro( << "Attempting to set a non-orthogonal matrix (after removing scaling)" );
    }

  typedef AdvancedMatrixOffsetTransformBase<TScalarType, 3> Baseclass;
  this->Baseclass::SetMatrix( matrix );
  this->PrecomputeJacobianOfSpatialJacobian();
}


// Set Parameters
template <class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>
::SetParameters( const ParametersType & parameters )
{

  itkDebugMacro( << "Setting parameters " << parameters );

  // Transfer the versor part
  
  AxisType axis;

  double norm = parameters[0]*parameters[0];
  axis[0] = parameters[0];
  norm += parameters[1]*parameters[1];
  axis[1] = parameters[1];
  norm += parameters[2]*parameters[2];
  axis[2] = parameters[2];
  if( norm > 0)
    {
    norm = vcl_sqrt(norm);
    }

  double epsilon = 1e-10;
  if(norm >= 1.0-epsilon)
    {
    axis = axis / (norm+epsilon*norm);
    }
  VersorType newVersor;
  newVersor.Set(axis);
  this->SetVarVersor( newVersor );
  m_Scale = parameters[6]; // must be set before calling ComputeMatrix();
  this->ComputeMatrix();

  itkDebugMacro( <<"Versor is now " << this->GetVersor() );
  
   
  // Transfer the translation part
  TranslationType newTranslation;
  newTranslation[0] = parameters[3];
  newTranslation[1] = parameters[4];
  newTranslation[2] = parameters[5];
  this->SetVarTranslation(newTranslation);
  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  itkDebugMacro(<<"After setting parameters ");
}


//
// Get Parameters
// 
// Parameters are ordered as:
//
// p[0:2] = right part of the versor (axis times vcl_sin(t/2))
// p[3:5} = translation components
// p[6:6} = scaling factor (isotropic)
//

template <class TScalarType>
const typename AdvancedSimilarity3DTransform<TScalarType>::ParametersType &
AdvancedSimilarity3DTransform<TScalarType>
::GetParameters( void ) const
{
  itkDebugMacro( << "Getting parameters ");

  this->m_Parameters[0] = this->GetVersor().GetX();
  this->m_Parameters[1] = this->GetVersor().GetY();
  this->m_Parameters[2] = this->GetVersor().GetZ();

  // Transfer the translation
  this->m_Parameters[3] = this->GetTranslation()[0];
  this->m_Parameters[4] = this->GetTranslation()[1];
  this->m_Parameters[5] = this->GetTranslation()[2];

  this->m_Parameters[6] = this->GetScale();

  itkDebugMacro(<<"After getting parameters " << this->m_Parameters );

  return this->m_Parameters;
}

// Set parameters
template<class TScalarType>
const typename AdvancedSimilarity3DTransform<TScalarType>::JacobianType &
AdvancedSimilarity3DTransform<TScalarType>::
GetJacobian( const InputPointType & p ) const
{
  /**
  typedef typename VersorType::ValueType  ValueType;
  

  // compute derivatives with respect to rotation
  const ValueType vx = this->GetVersor().GetX();
  const ValueType vy = this->GetVersor().GetY();
  const ValueType vz = this->GetVersor().GetZ();
  const ValueType vw = this->GetVersor().GetW();

  this->m_Jacobian.Fill(0.0);

  const InputVectorType pp = p - this->GetCenter();

  const double px = pp[0];
  const double py = pp[1];
  const double pz = pp[2];

  const double vxx = vx * vx;
  const double vyy = vy * vy;
  const double vzz = vz * vz;
  const double vww = vw * vw;

  const double vxy = vx * vy;
  const double vxz = vx * vz;
  const double vxw = vx * vw;

  const double vyz = vy * vz;
  const double vyw = vy * vw;

  const double vzw = vz * vw;


  // compute Jacobian with respect to quaternion parameters
  // SK: \todo: shoudn't these be multiplied by m_Scale?
  this->m_Jacobian[0][0] = 2.0 * (               (vyw+vxz)*py + (vzw-vxy)*pz)
                         / vw;
  this->m_Jacobian[1][0] = 2.0 * ((vyw-vxz)*px   -2*vxw   *py + (vxx-vww)*pz) 
                         / vw;
  this->m_Jacobian[2][0] = 2.0 * ((vzw+vxy)*px + (vww-vxx)*py   -2*vxw   *pz) 
                         / vw;

  this->m_Jacobian[0][1] = 2.0 * ( -2*vyw  *px + (vxw+vyz)*py + (vww-vyy)*pz) 
                         / vw;
  this->m_Jacobian[1][1] = 2.0 * ((vxw-vyz)*px                + (vzw+vxy)*pz) 
                         / vw;
  this->m_Jacobian[2][1] = 2.0 * ((vyy-vww)*px + (vzw-vxy)*py   -2*vyw   *pz) 
                         / vw;

  this->m_Jacobian[0][2] = 2.0 * ( -2*vzw  *px + (vzz-vww)*py + (vxw-vyz)*pz) 
                         / vw;
  this->m_Jacobian[1][2] = 2.0 * ((vww-vzz)*px   -2*vzw   *py + (vyw+vxz)*pz) 
                         / vw;
  this->m_Jacobian[2][2] = 2.0 * ((vxw+vyz)*px + (vyw-vxz)*py               ) 
                         / vw;
  */

  this->m_Jacobian.Fill(0.0);
  const InputVectorType pp = p - this->GetCenter();
  const JacobianOfSpatialJacobianType & jsj = this->m_JacobianOfSpatialJacobian;

  /** Compute dR/dmu * (p-c) */
  for(unsigned int dim=0; dim < SpaceDimension; dim++ ) 
  {
    const InputVectorType column = jsj[dim] * pp;
    for (unsigned int i=0; i < SpaceDimension; ++i)
    {
      this->m_Jacobian(i,dim) = column[i];
    }
  }

  // compute Jacobian with respect to the translation parameters
  this->m_Jacobian[0][3] = 1.0;
  this->m_Jacobian[1][4] = 1.0;
  this->m_Jacobian[2][5] = 1.0;

  // compute Jacobian with respect to the scale parameter
  const MatrixType & matrix = this->GetMatrix();

  const InputVectorType mpp = matrix * pp;

  this->m_Jacobian[0][6] = mpp[0] / m_Scale;
  this->m_Jacobian[1][6] = mpp[1] / m_Scale;
  this->m_Jacobian[2][6] = mpp[2] / m_Scale;

  return this->m_Jacobian;

}
 

// Set the scale factor
template <class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>
::ComputeMatrix()
{
  this->Superclass::ComputeMatrix();
  MatrixType newMatrix = this->GetMatrix();
  newMatrix *= m_Scale;
  this->SetVarMatrix(newMatrix);
  this->PrecomputeJacobianOfSpatialJacobian();
}


/** Compute the matrix */
template <class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>
::ComputeMatrixParameters( void )
{
  MatrixType matrix = this->GetMatrix();

  m_Scale = vnl_math_cuberoot( vnl_det( matrix.GetVnlMatrix() )  ) ;
  
  matrix /= m_Scale;

  VersorType v;
  v.Set( matrix );
  this->SetVarVersor( v );
  this->PrecomputeJacobianOfSpatialJacobian();

}

// Precompute Jacobian of Spatial Jacobian
template <class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>
::PrecomputeJacobianOfSpatialJacobian( void )
{
  if ( ParametersDimension < 7 )
  {
    return;
  }

  /** The Jacobian of spatial Jacobian remains constant, so is precomputed */
  JacobianOfSpatialJacobianType & jsj = this->m_JacobianOfSpatialJacobian;
  jsj.resize(ParametersDimension);

  typedef typename VersorType::ValueType  ValueType;

  // compute derivatives with respect to rotation
  const ValueType vx = this->GetVersor().GetX();
  const ValueType vy = this->GetVersor().GetY();
  const ValueType vz = this->GetVersor().GetZ();
  const ValueType vw = this->GetVersor().GetW();

  const double vxx = vx * vx;
  const double vyy = vy * vy;
  const double vzz = vz * vz;
  const double vww = vw * vw;

  const double vxy = vx * vy;
  const double vxz = vx * vz;
  const double vxw = vx * vw;

  const double vyz = vy * vz;
  const double vyw = vy * vw;

  const double vzw = vz * vw;

  jsj[0](0,0) = 0.0;     jsj[0](0,1) = vyw+vxz; jsj[0](0,2) = vzw-vxy;
  jsj[0](1,0) = vyw-vxz; jsj[0](1,1) =-2.0*vxw; jsj[0](1,2) = vxx-vww;
  jsj[0](2,0) = vzw+vxy; jsj[0](2,1) = vww-vxx; jsj[0](2,2) = -2.0*vxw ;
  jsj[0] *= ( this->m_Scale * 2.0 / vw );

  jsj[1](0,0) =-2.0*vyw; jsj[1](0,1) = vxw+vyz; jsj[1](0,2) = vww-vyy;
  jsj[1](1,0) = vxw-vyz; jsj[1](1,1) = 0.0;     jsj[1](1,2) = vzw+vxy;
  jsj[1](2,0) = vyy-vww; jsj[1](2,1) = vzw-vxy; jsj[1](2,2) =-2.0*vyw;
  jsj[1] *= ( this->m_Scale * 2.0 / vw );

  jsj[2](0,0) =-2.0*vzw; jsj[2](0,1) = vzz-vw;  jsj[2](0,2) = vxw-vyz;
  jsj[2](1,0) = vww-vzz; jsj[2](1,1) =-2.0*vzw; jsj[2](1,2) = vyw+vxz;
  jsj[2](2,0) = vxw+vyz; jsj[2](2,1) = vyw-vxz; jsj[2](2,2) = 0.0;
  jsj[2] *= ( this->m_Scale * 2.0 / vw );
  
  for ( unsigned int par = 3; par < 7; ++par )
  {
    jsj[par].Fill(0.0);
  }
  if ( vcl_abs( this->m_Scale ) > 0 )
  {
    jsj[6] = this->GetMatrix().GetVnlMatrix() / this->m_Scale;
  }
}

 
// Print self
template<class TScalarType>
void
AdvancedSimilarity3DTransform<TScalarType>::
PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Scale = " << m_Scale << std::endl;
}

} // namespace

#endif
