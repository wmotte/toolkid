/*======================================================================

  This file is part of the elastix software.

  Copyright (c) University Medical Center Utrecht. All rights reserved.
  See src/CopyrightElastix.txt or http://elastix.isi.uu.nl/legal.php for
  details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the above copyright notices for more information.

======================================================================*/

#ifndef __itkAdvancedCombinationTransform_hxx
#define __itkAdvancedCombinationTransform_hxx

#include "itkAdvancedCombinationTransform.h"


namespace itk
{
  
/**
 * ************************ Constructor *************************
 */

template <typename TScalarType, unsigned int NDimensions>
AdvancedCombinationTransform<TScalarType, NDimensions>
::AdvancedCombinationTransform() : Superclass(NDimensions,1)
{
  /** Initialize. */
  this->m_InitialTransform = 0;
  this->m_CurrentTransform = 0;

  /** Set composition by default.*/
  this->m_UseAddition = false;
  this->m_UseComposition = true;

  /** Set everything to have no current transform. */
  this->m_SelectedTransformPointFunction
    = &Self::TransformPointNoCurrentTransform;
  this->m_SelectedGetJacobianFunction
    = &Self::GetJacobianNoCurrentTransform;
  this->m_SelectedGetSparseJacobianFunction
    = &Self::GetJacobianNoCurrentTransform;
  this->m_SelectedGetSpatialJacobianFunction
    = &Self::GetSpatialJacobianNoCurrentTransform;
  this->m_SelectedGetSpatialHessianFunction
    = &Self::GetSpatialHessianNoCurrentTransform;
  this->m_SelectedGetJacobianOfSpatialJacobianFunction
    = &Self::GetJacobianOfSpatialJacobianNoCurrentTransform;
  this->m_SelectedGetJacobianOfSpatialJacobianFunction2
    = &Self::GetJacobianOfSpatialJacobianNoCurrentTransform;
  this->m_SelectedGetJacobianOfSpatialHessianFunction
    = &Self::GetJacobianOfSpatialHessianNoCurrentTransform;
  this->m_SelectedGetJacobianOfSpatialHessianFunction2
    = &Self::GetJacobianOfSpatialHessianNoCurrentTransform;

} // end Constructor
  

/**
 *
 * ***********************************************************
 * ***** Override functions to aid for combining transformations.
 *
 * ***********************************************************
 *
 */


/**
 * ***************** GetNumberOfParameters **************************
 */

template <typename TScalarType, unsigned int NDimensions>
unsigned int AdvancedCombinationTransform<TScalarType, NDimensions>
::GetNumberOfParameters( void ) const
{ 
  /** Return the number of parameters that completely define
   * the m_CurrentTransform.
   */
  if ( this->m_CurrentTransform.IsNotNull() )
  {
    return this->m_CurrentTransform->GetNumberOfParameters();
  }
  else
  {
    /** Throw an exception. */
    this->NoCurrentTransformSet();

    /** dummy return. */
    return this->m_Parameters.GetSize();
  }

} // end GetNumberOfParameters()


/**
 * ***************** GetNumberOfNonZeroJacobianIndices **************************
 */

template <typename TScalarType, unsigned int NDimensions>
unsigned long
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetNumberOfNonZeroJacobianIndices( void ) const
{ 
  /** Return the number of parameters that completely define
   * the m_CurrentTransform.
   */
  if ( this->m_CurrentTransform.IsNotNull() )
  {
    return this->m_CurrentTransform->GetNumberOfNonZeroJacobianIndices();
  }
  else
  {
    /** Throw an exception. */
    this->NoCurrentTransformSet();

    /** dummy return. */
    return this->m_Parameters.GetSize();
  }

} // end GetNumberOfNonZeroJacobianIndices()


/**
 * ***************** IsLinear **************************
 */

template <typename TScalarType, unsigned int NDimensions>
bool AdvancedCombinationTransform<TScalarType, NDimensions>
::IsLinear( void ) const
{ 
  bool currentLinear = true;
  if ( this->m_CurrentTransform.IsNotNull() )
  {
    currentLinear = this->m_CurrentTransform->IsLinear();
  }

  bool initialLinear = true;
  if ( this->m_InitialTransform.IsNotNull() )
  {
    initialLinear = this->m_InitialTransform->IsLinear();
  }

  return ( currentLinear && initialLinear );

} // end IsLinear()


/**
 * ***************** GetParameters **************************
 */

template <typename TScalarType, unsigned int NDimensions>
const typename AdvancedCombinationTransform<TScalarType, NDimensions>::ParametersType &
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetParameters( void ) const
{ 
  /** Return the parameters that completely define the m_CurrentTransform. */
  if ( this->m_CurrentTransform.IsNotNull() )
  {
    return this->m_CurrentTransform->GetParameters();
  }
  else
  {
    /** Throw an exception. */
    this->NoCurrentTransformSet();
    
    /** dummy return. */
    return this->m_Parameters;
  }

} // end GetParameters()


/**
 * ***************** SetParameters **************************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::SetParameters( const ParametersType & param )
{ 
  /** Set the parameters in the m_CurrentTransform. */
  if ( this->m_CurrentTransform.IsNotNull() )
  {
    this->Modified();
    this->m_CurrentTransform->SetParameters( param );
  }
  else
  {
    /** Throw an exception. */
    this->NoCurrentTransformSet();
  }

} // end SetParameters()


/**
 * ***************** SetParametersByValue **************************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::SetParametersByValue( const ParametersType & param )
{ 
  /** Set the parameters in the m_CurrentTransfom. */
  if ( this->m_CurrentTransform.IsNotNull() )
  {
    this->Modified();
    this->m_CurrentTransform->SetParametersByValue( param );
  }
  else
  {
    /** Throw an exception. */
    this->NoCurrentTransformSet();
  }

} // end SetParametersByValue()


/**
 * ***************** GetInverse **************************
 */

template <typename TScalarType, unsigned int NDimensions>
bool AdvancedCombinationTransform<TScalarType, NDimensions>
::GetInverse( Self * inverse ) const
{ 
  if( !inverse )
  {
    /** Inverse transformation cannot be returned into nothingness. */
    return false;
  }
  else if ( this->m_CurrentTransform.IsNull() )
  { 
    /** No current transform has been set. Throw an exception. */
    this->NoCurrentTransformSet();
    return false;
  }
  else if ( this->m_InitialTransform.IsNull() )
  {
    /** No Initial transform, so call the CurrentTransform's implementation. */
    return this->m_CurrentTransform->GetInverse( inverse );
  }
  else if ( this->m_UseAddition )
  {
    /** No generic expression exists for the inverse of (T0+T1)(x). */
    return false;
  }
  else // UseComposition
  {
    /** The initial transform and the current transform have been set
     * and UseComposition is set to true.
     * The inverse transform IT is defined by:
     *  IT ( T1(T0(x) ) = x
     * So:
     *  IT(y) = T0^{-1} ( T1^{-1} (y) ),
     * which is of course only defined when the inverses of both
     * the initial and the current transforms are defined.
     */

    /** Try create the inverse of the initial transform. */
    InitialTransformPointer inverseT0 = InitialTransformType::New();
    bool T0invertable = this->m_InitialTransform->GetInverse( inverseT0 );

    if ( T0invertable )
    {
      /** Try to create the inverse of the current transform. */
      CurrentTransformPointer inverseT1 = CurrentTransformType::New();
      bool T1invertable = this->m_CurrentTransform->GetInverse( inverseT1 );

      if ( T1invertable )
      {
        /** The transform can be inverted! */
        inverse->SetUseComposition( true );
        inverse->SetInitialTransform( inverseT1 );
        inverse->SetCurrentTransform( inverseT0 );
        return true;
      }
      else
      {
        /** The initial transform is invertible, but the current one not. */
        return false;
      }
    }
    else
    {
      /** The initial transform is not invertible. */
      return false;
    }

  } // end else: UseComposition

} // end GetInverse()


/**
 * ***************** GetHasNonZeroSpatialHessian **************************
 */

template <typename TScalarType, unsigned int NDimensions>
bool
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetHasNonZeroSpatialHessian( void ) const
{
  /** Set the parameters in the m_CurrentTransfom. */
  if ( this->m_CurrentTransform.IsNull() )
  {
    /** No current transform has been set. Throw an exception. */
    this->NoCurrentTransformSet();
    return false;
  }
  else if ( this->m_InitialTransform.IsNull() )
  {
    /** No Initial transform, so call the CurrentTransform's implementation. */
    return this->m_CurrentTransform->GetHasNonZeroSpatialHessian();
  }
  else
  {
    bool dummy = this->m_InitialTransform->GetHasNonZeroSpatialHessian()
      || this->m_CurrentTransform->GetHasNonZeroSpatialHessian();
    return dummy;
  }

} // end GetHasNonZeroSpatialHessian()


/**
 * ***************** HasNonZeroJacobianOfSpatialHessian **************************
 */

template <typename TScalarType, unsigned int NDimensions>
bool
AdvancedCombinationTransform<TScalarType, NDimensions>
::HasNonZeroJacobianOfSpatialHessian( void ) const
{
  /** Set the parameters in the m_CurrentTransfom. */
  if ( this->m_CurrentTransform.IsNull() )
  {
    /** No current transform has been set. Throw an exception. */
    this->NoCurrentTransformSet();
    return false;
  }
  else if ( this->m_InitialTransform.IsNull() )
  {
    /** No Initial transform, so call the CurrentTransform's implementation. */
    return this->m_CurrentTransform->GetHasNonZeroJacobianOfSpatialHessian();
  }
  else
  {
    bool dummy = this->m_InitialTransform->GetHasNonZeroJacobianOfSpatialHessian()
      || this->m_CurrentTransform->GetHasNonZeroJacobianOfSpatialHessian();
    return dummy;
  }

} // end HasNonZeroJacobianOfSpatialHessian()


/**
 *
 * ***********************************************************
 * ***** Functions to set the transformations and choose the
 * ***** combination method.
 *
 * ***********************************************************
 *
 */


/**
 * ******************* SetInitialTransform **********************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::SetInitialTransform( const InitialTransformType * _arg )
{
  /** Set the the initial transform and call the UpdateCombinationMethod. */
  if ( this->m_InitialTransform != _arg )
  {
    this->m_InitialTransform = _arg;
    this->Modified();
    this->UpdateCombinationMethod();
  }

} // end SetInitialTransform()


/**
 * ******************* SetCurrentTransform **********************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::SetCurrentTransform( CurrentTransformType * _arg )
{
  /** Set the the current transform and call the UpdateCombinationMethod. */
  if ( this->m_CurrentTransform != _arg )
  {
    this->m_CurrentTransform = _arg;
    this->Modified();
    this->UpdateCombinationMethod();
  }

} // end SetCurrentTransform()


/**
 * ********************** SetUseAddition **********************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::SetUseAddition( bool _arg )
{
  /** Set the UseAddition and UseComposition bools and call the UpdateCombinationMethod. */
  if ( this->m_UseAddition != _arg )
  {
    this->m_UseAddition = _arg;
    this->m_UseComposition = !_arg;
    this->Modified();
    this->UpdateCombinationMethod();
  }

} // end SetUseAddition()


/**
 * ********************** SetUseComposition *******************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::SetUseComposition( bool _arg )
{
  /** Set the UseAddition and UseComposition bools and call the UpdateCombinationMethod. */
  if ( this->m_UseComposition != _arg )
  {
    this->m_UseComposition = _arg;
    this->m_UseAddition = !_arg;
    this->Modified();
    this->UpdateCombinationMethod();
  }

} // end SetUseComposition()


/**
 * ****************** UpdateCombinationMethod ********************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::UpdateCombinationMethod( void )
{
  /** Update the m_SelectedTransformPointFunction and
   * the m_SelectedGetJacobianFunction
   */
  if ( this->m_CurrentTransform.IsNull() )
  {
    this->m_SelectedTransformPointFunction
      = &Self::TransformPointNoCurrentTransform;
    this->m_SelectedGetJacobianFunction
      = &Self::GetJacobianNoCurrentTransform;
    this->m_SelectedGetSparseJacobianFunction
      = &Self::GetJacobianNoCurrentTransform;
    this->m_SelectedGetSpatialJacobianFunction
      = &Self::GetSpatialJacobianNoCurrentTransform;
    this->m_SelectedGetSpatialHessianFunction
      = &Self::GetSpatialHessianNoCurrentTransform;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction
      = &Self::GetJacobianOfSpatialJacobianNoCurrentTransform;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction2
      = &Self::GetJacobianOfSpatialJacobianNoCurrentTransform;
    this->m_SelectedGetJacobianOfSpatialHessianFunction
      = &Self::GetJacobianOfSpatialHessianNoCurrentTransform;
    this->m_SelectedGetJacobianOfSpatialHessianFunction2
      = &Self::GetJacobianOfSpatialHessianNoCurrentTransform;
  }
  else if ( this->m_InitialTransform.IsNull() )
  {
    this->m_SelectedTransformPointFunction
      = &Self::TransformPointNoInitialTransform;
    this->m_SelectedGetJacobianFunction
      = &Self::GetJacobianNoInitialTransform;
    this->m_SelectedGetSparseJacobianFunction
      = &Self::GetJacobianNoInitialTransform;
    this->m_SelectedGetSpatialJacobianFunction
      = &Self::GetSpatialJacobianNoInitialTransform;
    this->m_SelectedGetSpatialHessianFunction
      = &Self::GetSpatialHessianNoInitialTransform;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction
      = &Self::GetJacobianOfSpatialJacobianNoInitialTransform;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction2
      = &Self::GetJacobianOfSpatialJacobianNoInitialTransform;
    this->m_SelectedGetJacobianOfSpatialHessianFunction
      = &Self::GetJacobianOfSpatialHessianNoInitialTransform;
    this->m_SelectedGetJacobianOfSpatialHessianFunction2
      = &Self::GetJacobianOfSpatialHessianNoInitialTransform;
  }
  else if ( this->m_UseAddition )
  {
    this->m_SelectedTransformPointFunction
      = &Self::TransformPointUseAddition;
    this->m_SelectedGetJacobianFunction
      = &Self::GetJacobianUseAddition;
    this->m_SelectedGetSparseJacobianFunction
      = &Self::GetJacobianUseAddition;
    this->m_SelectedGetSpatialJacobianFunction
      = &Self::GetSpatialJacobianUseAddition;
    this->m_SelectedGetSpatialHessianFunction
      = &Self::GetSpatialHessianUseAddition;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction
      = &Self::GetJacobianOfSpatialJacobianUseAddition;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction2
      = &Self::GetJacobianOfSpatialJacobianUseAddition;
    this->m_SelectedGetJacobianOfSpatialHessianFunction
      = &Self::GetJacobianOfSpatialHessianUseAddition;
    this->m_SelectedGetJacobianOfSpatialHessianFunction2
      = &Self::GetJacobianOfSpatialHessianUseAddition;
  }
  else
  {
    this->m_SelectedTransformPointFunction
      = &Self::TransformPointUseComposition;
    this->m_SelectedGetJacobianFunction
      = &Self::GetJacobianUseComposition;
    this->m_SelectedGetSparseJacobianFunction
      = &Self::GetJacobianUseComposition;
    this->m_SelectedGetSpatialJacobianFunction
      = &Self::GetSpatialJacobianUseComposition;
    this->m_SelectedGetSpatialHessianFunction
      = &Self::GetSpatialHessianUseComposition;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction
      = &Self::GetJacobianOfSpatialJacobianUseComposition;
    this->m_SelectedGetJacobianOfSpatialJacobianFunction2
      = &Self::GetJacobianOfSpatialJacobianUseComposition;
    this->m_SelectedGetJacobianOfSpatialHessianFunction
      = &Self::GetJacobianOfSpatialHessianUseComposition;
    this->m_SelectedGetJacobianOfSpatialHessianFunction2
      = &Self::GetJacobianOfSpatialHessianUseComposition;
  }

} // end UpdateCombinationMethod()


/**
 * ************* NoCurrentTransformSet **********************
 */

template <typename TScalarType, unsigned int NDimensions>
void AdvancedCombinationTransform<TScalarType, NDimensions>
::NoCurrentTransformSet( void ) const throw (ExceptionObject)
{
  itkExceptionMacro( << "No current transform set in the AdvancedCombinationTransform" );

} // end NoCurrentTransformSet()


/**
 *
 * ***********************************************************
 * ***** Functions that implement the:
 * ***** - TransformPoint()
 * ***** - GetJacobian()
 * ***** - GetSpatialJacobian()
 * ***** - GetSpatialHessian()
 * ***** - GetJacobianOfSpatialJacobian()
 * ***** - GetJacobianOfSpatialHessian()
 * ***** for the four possible cases:
 * ***** - no initial transform: this is the same as using only one transform
 * ***** - no current transform: error, it should be set
 * ***** - use addition to combine transformations
 * ***** - use composition to combine transformations
 *
 * ***********************************************************
 *
 */


/**
 * ************* TransformPointUseAddition **********************
 */

template <typename TScalarType, unsigned int NDimensions>
typename AdvancedCombinationTransform<TScalarType, NDimensions>::OutputPointType
AdvancedCombinationTransform<TScalarType, NDimensions>
::TransformPointUseAddition( const InputPointType & point ) const
{       
  /** The Initial transform. */
  OutputPointType out0
    = this->m_InitialTransform->TransformPoint( point );

  /** The Current transform. */
    OutputPointType out
      = this->m_CurrentTransform->TransformPoint( point );

  /** Add them. */
  for ( unsigned int i=0; i < SpaceDimension; i++ )
  {
    out[ i ] += ( out0[ i ] - point[ i ] );
  }

  return out;

} // end TransformPointUseAddition()


/**
 * **************** TransformPointUseComposition *************
 */

template <typename TScalarType, unsigned int NDimensions>
typename AdvancedCombinationTransform<TScalarType, NDimensions>::OutputPointType
AdvancedCombinationTransform<TScalarType, NDimensions>
::TransformPointUseComposition( const InputPointType & point ) const
{
  return this->m_CurrentTransform->TransformPoint(
    this->m_InitialTransform->TransformPoint( point ) );

} // end TransformPointUseComposition()


/**
 * **************** TransformPointNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
typename AdvancedCombinationTransform<TScalarType, NDimensions>::OutputPointType
AdvancedCombinationTransform<TScalarType, NDimensions>
::TransformPointNoInitialTransform( const InputPointType & point ) const
{
  return this->m_CurrentTransform->TransformPoint( point );

} // end TransformPointNoInitialTransform()


/**
 * ******** TransformPointNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
typename AdvancedCombinationTransform<TScalarType, NDimensions>::OutputPointType
AdvancedCombinationTransform<TScalarType, NDimensions>
::TransformPointNoCurrentTransform( const InputPointType & point ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();
  
  /** dummy return. */
  return point;

} // end TransformPointNoCurrentTransform()


/**
 * ************* GetJacobianUseAddition ***************************
 */

template <typename TScalarType, unsigned int NDimensions>
const typename AdvancedCombinationTransform<TScalarType, NDimensions>::JacobianType &
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianUseAddition( const InputPointType & point ) const
{
  return this->m_CurrentTransform->GetJacobian( point );

} // end GetJacobianUseAddition()


/**
 * **************** GetJacobianUseComposition *************
 */

template <typename TScalarType, unsigned int NDimensions>
const typename AdvancedCombinationTransform<TScalarType, NDimensions>::JacobianType &
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianUseComposition( const InputPointType & point ) const
{
  return this->m_CurrentTransform->GetJacobian(
    this->m_InitialTransform->TransformPoint( point ) );

} // end GetJacobianUseComposition()


/**
 * **************** GetJacobianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
const typename AdvancedCombinationTransform<TScalarType, NDimensions>::JacobianType &
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianNoInitialTransform( const InputPointType & point ) const
{
  return this->m_CurrentTransform->GetJacobian( point );

} // end GetJacobianNoInitialTransform()


/**
 * ******** GetJacobianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
const typename AdvancedCombinationTransform<TScalarType, NDimensions>::JacobianType &
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianNoCurrentTransform( const InputPointType & point ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

  /** dummy return. */
  return this->m_Jacobian;

} // end GetJacobianNoCurrentTransform()


/**
 * ************* GetJacobianUseAddition ***************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianUseAddition(
  const InputPointType & ipp,
  JacobianType & j,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobian( ipp, j, nonZeroJacobianIndices );

} // end GetJacobianUseAddition()


/**
 * **************** GetJacobianUseComposition *************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianUseComposition(
  const InputPointType & ipp,
  JacobianType & j,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobian(
    this->m_InitialTransform->TransformPoint( ipp ),
    j, nonZeroJacobianIndices );

} // end GetJacobianUseComposition()


/**
 * **************** GetJacobianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianNoInitialTransform(
  const InputPointType & ipp,
  JacobianType & j,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobian( ipp, j, nonZeroJacobianIndices );

} // end GetJacobianNoInitialTransform()


/**
 * ******** GetJacobianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianNoCurrentTransform(
  const InputPointType & ipp,
  JacobianType & j,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetJacobianNoCurrentTransform()


/**
 * ************* GetSpatialJacobianUseAddition ***************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialJacobianUseAddition(
  const InputPointType & ipp,
  SpatialJacobianType & sj ) const
{
  SpatialJacobianType sj0, sj1, identity;
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_CurrentTransform->GetSpatialJacobian( ipp, sj1 );
  identity.SetIdentity();
  sj = sj0 + sj1 - identity;

} // end GetSpatialJacobianUseAddition()


/**
 * **************** GetSpatialJacobianUseComposition *************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialJacobianUseComposition(
  const InputPointType & ipp,
  SpatialJacobianType & sj ) const
{
  SpatialJacobianType sj0, sj1;
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_CurrentTransform->GetSpatialJacobian(
    this->m_InitialTransform->TransformPoint( ipp ), sj1 );

  sj = sj1 * sj0;    

} // end GetSpatialJacobianUseComposition()


/**
 * **************** GetSpatialJacobianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialJacobianNoInitialTransform(
  const InputPointType & ipp,
  SpatialJacobianType & sj ) const
{
  this->m_CurrentTransform->GetSpatialJacobian( ipp, sj );

} // end GetSpatialJacobianNoInitialTransform()


/**
 * ******** GetSpatialJacobianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialJacobianNoCurrentTransform(
  const InputPointType & ipp,
  SpatialJacobianType & sj ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetSpatialJacobianNoCurrentTransform()


/**
 * ******** GetSpatialHessianUseAddition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialHessianUseAddition(
  const InputPointType & ipp,
  SpatialHessianType & sh ) const
{
  SpatialHessianType sh0, sh1;
  this->m_InitialTransform->GetSpatialHessian( ipp, sh0 );
  this->m_CurrentTransform->GetSpatialHessian( ipp, sh1 );

  for ( unsigned int i = 0; i < SpaceDimension; ++i )
  {
    sh[ i ] = sh0[ i ] + sh1[ i ];
  }

} // end GetSpatialHessianUseAddition()


/**
 * ******** GetSpatialHessianUseComposition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialHessianUseComposition(
  const InputPointType & ipp,
  SpatialHessianType & sh ) const
{
  /** Create intermediary variables for the internal transforms. */
  SpatialJacobianType sj0, sj1;
  SpatialHessianType sh0, sh1;

  /** Transform the input point. */
  // \todo this has already been computed and it is expensive.
  InputPointType transformedPoint
    = this->m_InitialTransform->TransformPoint( ipp );

  /** Compute the (Jacobian of the) spatial Jacobian / Hessian of the
   * internal transforms.
   */
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_CurrentTransform->GetSpatialJacobian( transformedPoint, sj1 );
  this->m_InitialTransform->GetSpatialHessian( ipp, sh0 );
  this->m_CurrentTransform->GetSpatialHessian( transformedPoint, sh1 );

  typename SpatialJacobianType::InternalMatrixType sj0tvnl = sj0.GetTranspose();
  SpatialJacobianType sj0t( sj0tvnl );

  /** Combine them in one overall spatial Hessian. */
  for ( unsigned int dim = 0; dim < SpaceDimension; ++dim )
  {
    sh[dim] = sj0t * ( sh1[dim] * sj0 );

    for ( unsigned int p = 0; p < SpaceDimension; ++p )
    {
      sh[dim] += ( sh0[p] * sj1(dim,p) ); 
    }
  }   

} // end GetSpatialHessianUseComposition()


/**
 * ******** GetSpatialHessianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialHessianNoInitialTransform(
  const InputPointType & ipp,
  SpatialHessianType & sh ) const
{
  this->m_CurrentTransform->GetSpatialHessian( ipp, sh );

} // end GetSpatialHessianNoInitialTransform()


/**
 * ******** GetSpatialHessianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialHessianNoCurrentTransform(
  const InputPointType & ipp,
  SpatialHessianType & sh ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetSpatialHessianNoCurrentTransform()


/**
 * ******** GetJacobianOfSpatialJacobianUseAddition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianUseAddition(
  const InputPointType & ipp,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    ipp, jsj, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialJacobianUseAddition()


/**
 * ******** GetJacobianOfSpatialJacobianUseAddition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianUseAddition(
  const InputPointType & ipp,
  SpatialJacobianType & sj,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    ipp, sj, jsj, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialJacobianUseAddition()


/**
 * ******** GetJacobianOfSpatialJacobianUseComposition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianUseComposition(
  const InputPointType & ipp,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  SpatialJacobianType sj0;
  JacobianOfSpatialJacobianType jsj1;
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    this->m_InitialTransform->TransformPoint( ipp ),
    jsj1, nonZeroJacobianIndices );

  jsj.resize( nonZeroJacobianIndices.size() );
  for ( unsigned int mu = 0; mu < nonZeroJacobianIndices.size(); ++mu )
  {
    jsj[ mu ] = jsj1[ mu ] * sj0;    
  }

} // end GetJacobianOfSpatialJacobianUseComposition()


/**
 * ******** GetJacobianOfSpatialJacobianUseComposition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianUseComposition(
  const InputPointType & ipp,
  SpatialJacobianType & sj,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  SpatialJacobianType sj0, sj1;
  JacobianOfSpatialJacobianType jsj1;
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    this->m_InitialTransform->TransformPoint( ipp ),
    sj1, jsj1, nonZeroJacobianIndices );
  
  sj = sj1 * sj0;    
  jsj.resize( nonZeroJacobianIndices.size() );
  for ( unsigned int mu = 0; mu < nonZeroJacobianIndices.size(); ++mu )
  {    
    jsj[ mu ] = jsj1[ mu ] * sj0;    
  }

} // end GetJacobianOfSpatialJacobianUseComposition()


/**
 * ******** GetJacobianOfSpatialJacobianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianNoInitialTransform(
  const InputPointType & ipp,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    ipp, jsj, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialJacobianNoInitialTransform()


/**
 * ******** GetJacobianOfSpatialJacobianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianNoInitialTransform(
  const InputPointType & ipp,
  SpatialJacobianType & sj,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    ipp, sj, jsj, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialJacobianNoInitialTransform()


/**
 * ******** GetJacobianOfSpatialJacobianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianNoCurrentTransform(
  const InputPointType & ipp,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetJacobianOfSpatialJacobianNoCurrentTransform()


/**
 * ******** GetJacobianOfSpatialJacobianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobianNoCurrentTransform(
  const InputPointType & ipp,
  SpatialJacobianType & sj,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetJacobianOfSpatialJacobianNoCurrentTransform()


/**
 * ******** GetJacobianOfSpatialHessianUseAddition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianUseAddition(
  const InputPointType & ipp,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialHessian(
    ipp, jsh, nonZeroJacobianIndices );
  
} // end GetJacobianOfSpatialHessianUseAddition()


/**
 * ******** GetJacobianOfSpatialHessianUseAddition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianUseAddition(
  const InputPointType & ipp,
  SpatialHessianType & sh,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialHessian(
    ipp, sh, jsh, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialHessianUseAddition()


/**
 * ******** GetJacobianOfSpatialHessianUseComposition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianUseComposition(
  const InputPointType & ipp,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Create intermediary variables for the internal transforms. */
  SpatialJacobianType sj0;
  SpatialHessianType sh0;
  JacobianOfSpatialJacobianType jsj1;
  JacobianOfSpatialHessianType jsh1;

  /** Transform the input point. */
  // \todo: this has already been computed and it is expensive.
  InputPointType transformedPoint
    = this->m_InitialTransform->TransformPoint( ipp );

  /** Compute the (Jacobian of the) spatial Jacobian / Hessian of the
   * internal transforms. */  
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_InitialTransform->GetSpatialHessian( ipp, sh0 );

  /** Assume/demand that GetJacobianOfSpatialJacobian returns
   * the same nonZeroJacobianIndices as the GetJacobianOfSpatialHessian. */
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    transformedPoint, jsj1, nonZeroJacobianIndices );
  this->m_CurrentTransform->GetJacobianOfSpatialHessian(
    transformedPoint, jsh1, nonZeroJacobianIndices );

  typename SpatialJacobianType::InternalMatrixType sj0tvnl = sj0.GetTranspose();
  SpatialJacobianType sj0t( sj0tvnl );
  
  jsh.resize( nonZeroJacobianIndices.size() );  

  /** Combine them in one overall Jacobian of spatial Hessian. */
  for ( unsigned int mu = 0; mu < nonZeroJacobianIndices.size(); ++mu )
  {
    for ( unsigned int dim = 0; dim < SpaceDimension; ++dim )
    {
      jsh[mu][dim] = sj0t * ( jsh1[mu][dim] * sj0 );

      for ( unsigned int p = 0; p < SpaceDimension; ++p )
      {
        jsh[mu][dim] += ( sh0[p] * jsj1[mu](dim,p) );
      }
    }
  }

} // end GetJacobianOfSpatialHessianUseComposition()


/**
 * ******** GetJacobianOfSpatialHessianUseComposition ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianUseComposition(
  const InputPointType & ipp,
  SpatialHessianType & sh,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Create intermediary variables for the internal transforms. */
  SpatialJacobianType sj0, sj1;
  SpatialHessianType sh0, sh1;
  JacobianOfSpatialJacobianType jsj1;
  JacobianOfSpatialHessianType jsh1;

  /** Transform the input point. */
  // \todo this has already been computed and it is expensive.
  InputPointType transformedPoint
    = this->m_InitialTransform->TransformPoint( ipp );

  /** Compute the (Jacobian of the) spatial Jacobian / Hessian of the
   * internal transforms. */
  this->m_InitialTransform->GetSpatialJacobian( ipp, sj0 );
  this->m_InitialTransform->GetSpatialHessian( ipp, sh0 );

    /** Assume/demand that GetJacobianOfSpatialJacobian returns
   * the same nonZeroJacobianIndices as the GetJacobianOfSpatialHessian. */
  this->m_CurrentTransform->GetJacobianOfSpatialJacobian(
    transformedPoint, sj1, jsj1, nonZeroJacobianIndices );
  this->m_CurrentTransform->GetJacobianOfSpatialHessian(
    transformedPoint, sh1, jsh1, nonZeroJacobianIndices );

  typename SpatialJacobianType::InternalMatrixType sj0tvnl = sj0.GetTranspose();
  SpatialJacobianType sj0t( sj0tvnl );  
  jsh.resize( nonZeroJacobianIndices.size() );  

  /** Combine them in one overall Jacobian of spatial Hessian. */
  for ( unsigned int mu = 0; mu < nonZeroJacobianIndices.size(); ++mu )
  {
    for ( unsigned int dim = 0; dim < SpaceDimension; ++dim )
    {
      jsh[mu][dim] = sj0t * ( jsh1[mu][dim] * sj0 );

      for ( unsigned int p = 0; p < SpaceDimension; ++p )
      {
        jsh[mu][dim] += ( sh0[p] * jsj1[mu](dim,p) );
      }
    }
  }

   /** Combine them in one overall spatial Hessian. */
  for ( unsigned int dim = 0; dim < SpaceDimension; ++dim )
  {
    sh[dim] = sj0t * ( sh1[dim] * sj0 );

    for ( unsigned int p = 0; p < SpaceDimension; ++p )
    {
      sh[dim] += ( sh0[p] * sj1(dim,p) );
    }
  }   

} // end GetJacobianOfSpatialHessianUseComposition()


/**
 * ******** GetJacobianOfSpatialHessianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianNoInitialTransform(
  const InputPointType & ipp,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialHessian(
    ipp, jsh, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialHessianNoInitialTransform()


/**
 * ******** GetJacobianOfSpatialHessianNoInitialTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianNoInitialTransform(
  const InputPointType & ipp,
  SpatialHessianType & sh,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  this->m_CurrentTransform->GetJacobianOfSpatialHessian(
    ipp, sh, jsh, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialHessianNoInitialTransform()


/**
 * ******** GetJacobianOfSpatialHessianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianNoCurrentTransform(
  const InputPointType & ipp,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetJacobianOfSpatialHessianNoCurrentTransform()


/**
 * ******** GetJacobianOfSpatialHessianNoCurrentTransform ******************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessianNoCurrentTransform(
  const InputPointType & ipp,
  SpatialHessianType & sh,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Throw an exception. */
  this->NoCurrentTransformSet();

} // end GetJacobianOfSpatialHessianNoCurrentTransform()


/**
 *
 * ***********************************************************
 * ***** Functions that point to the selected implementation.
 *
 * ***********************************************************
 *
 */


/**
 * ****************** TransformPoint ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
typename AdvancedCombinationTransform<TScalarType, NDimensions>::OutputPointType
AdvancedCombinationTransform<TScalarType, NDimensions>
::TransformPoint( const InputPointType & point ) const
{ 
  /** Call the selected TransformPointFunction. */
  return ((*this).*m_SelectedTransformPointFunction)( point );

} // end TransformPoint()


/**
 * ****************** GetJacobian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
const typename AdvancedCombinationTransform<TScalarType, NDimensions>::JacobianType &
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobian( const InputPointType & point ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetJacobianFunction)( point );

} // end GetJacobian()


/**
 * ****************** GetJacobian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobian(
  const InputPointType & ipp,
  JacobianType & j,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetSparseJacobianFunction)(
    ipp, j, nonZeroJacobianIndices );

} // end GetJacobian()


/**
 * ****************** GetSpatialJacobian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialJacobian(
  const InputPointType & ipp,
  SpatialJacobianType & sj ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetSpatialJacobianFunction)(
    ipp, sj );

} // end GetSpatialJacobian()


/**
 * ****************** GetSpatialHessian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetSpatialHessian(
  const InputPointType & ipp,
  SpatialHessianType & sh ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetSpatialHessianFunction)(
    ipp, sh );

} // end GetSpatialHessian()


/**
 * ****************** GetJacobianOfSpatialJacobian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobian(
  const InputPointType & ipp,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetJacobianOfSpatialJacobianFunction)(
    ipp, jsj, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialJacobian()


/**
 * ****************** GetJacobianOfSpatialJacobian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialJacobian(
  const InputPointType & ipp,
  SpatialJacobianType & sj,
  JacobianOfSpatialJacobianType & jsj,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetJacobianOfSpatialJacobianFunction2)(
    ipp, sj, jsj, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialJacobian()


/**
 * ****************** GetJacobianOfSpatialHessian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessian(
  const InputPointType & ipp,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetJacobianOfSpatialHessianFunction)(
    ipp, jsh, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialHessian()


/**
 * ****************** GetJacobianOfSpatialHessian ****************************
 */

template <typename TScalarType, unsigned int NDimensions>
void
AdvancedCombinationTransform<TScalarType, NDimensions>
::GetJacobianOfSpatialHessian(
  const InputPointType & ipp,
  SpatialHessianType & sh,
  JacobianOfSpatialHessianType & jsh,
  NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const
{
  /** Call the selected GetJacobian. */
  return ((*this).*m_SelectedGetJacobianOfSpatialHessianFunction2)(
    ipp, sh, jsh, nonZeroJacobianIndices );

} // end GetJacobianOfSpatialHessian()


} // end namespace itk


#endif // end #ifndef __itkAdvancedCombinationTransform_hxx

