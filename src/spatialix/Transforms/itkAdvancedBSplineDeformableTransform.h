/*======================================================================

This file is part of the elastix software.

Copyright (c) University Medical Center Utrecht. All rights reserved.
See src/CopyrightElastix.txt or http://elastix.isi.uu.nl/legal.php for
details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the above copyright notices for more information.

======================================================================*/

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdvancedBSplineDeformableTransform.h,v $
  Language:  C++
  Date:      $Date: 2008-04-11 16:28:11 $
  Version:   $Revision: 1.38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkAdvancedBSplineDeformableTransform_h
#define __itkAdvancedBSplineDeformableTransform_h

#include "itkAdvancedTransform.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkBSplineInterpolationWeightFunction2.h"
#include "itkBSplineInterpolationDerivativeWeightFunction.h"
#include "itkBSplineInterpolationSecondOrderDerivativeWeightFunction.h"


namespace itk
{

/** \class AdvancedBSplineDeformableTransform
 * \brief Deformable transform using a BSpline representation
 *
 * This class encapsulates a deformable transform of points from one 
 * N-dimensional one space to another N-dimensional space.
 * The deformation field is modeled using B-splines. 
 * A deformation is defined on a sparse regular grid of control points
 * \f$ \vec{\lambda}_j \f$ and is varied by defining a deformation 
 * \f$ \vec{g}(\vec{\lambda}_j) \f$ of each control point. 
 * The deformation \f$ D(\vec{x}) \f$ at any point \f$ \vec{x} \f$
 * is obtained by using a B-spline interpolation kernel.
 *
 * The deformation field grid is defined by a user specified GridRegion, 
 * GridSpacing and GridOrigin. Each grid/control point has associated with it
 * N deformation coefficients \f$ \vec{\delta}_j \f$, representing the N 
 * directional components of the deformation. Deformation outside the grid
 * plus support region for the BSpline interpolation is assumed to be zero.
 *
 * Additionally, the user can specified an addition bulk transform \f$ B \f$
 * such that the transformed point is given by:
 * \f[ \vec{y} = B(\vec{x}) + D(\vec{x}) \f]
 *
 * The parameters for this transform is N x N-D grid of spline coefficients.
 * The user specifies the parameters as one flat array: each N-D grid
 * is represented by an array in the same way an N-D image is represented
 * in the buffer; the N arrays are then concatentated together on form
 * a single array.
 *
 * For efficiency, this transform does not make a copy of the parameters.
 * It only keeps a pointer to the input parameters and assumes that the memory
 * is managed by the caller.
 *
 * The following illustrates the typical usage of this class:
 * \verbatim
 * typedef AdvancedBSplineDeformableTransform<double,2,3> TransformType;
 * TransformType::Pointer transform = TransformType::New();
 *
 * transform->SetGridRegion( region );
 * transform->SetGridSpacing( spacing );
 * transform->SetGridOrigin( origin );
 *
 * // NB: the region must be set first before setting the parameters
 *
 * TransformType::ParametersType parameters( 
 *                                       transform->GetNumberOfParameters() );
 *
 * // Fill the parameters with values
 *
 * transform->SetParameters( parameters )
 *
 * outputPoint = transform->TransformPoint( inputPoint );
 *
 * \endverbatim
 *
 * An alternative way to set the B-spline coefficients is via array of
 * images. The grid region, spacing and origin information is taken
 * directly from the first image. It is assumed that the subsequent images
 * are the same buffered region. The following illustrates the API:
 * \verbatim
 *
 * TransformType::ImageConstPointer images[2];
 *
 * // Fill the images up with values
 *
 * transform->SetCoefficientImages( images ); 
 * outputPoint = transform->TransformPoint( inputPoint );
 *
 * \endverbatim
 *
 * Warning: use either the SetParameters() or SetCoefficientImage()
 * API. Mixing the two modes may results in unexpected results.
 *
 * The class is templated coordinate representation type (float or double),
 * the space dimension and the spline order.
 *
 * \ingroup Transforms
 */
template <
    class TScalarType = double,          // Data type for scalars
    unsigned int NDimensions = 3,        // Number of dimensions
    unsigned int VSplineOrder = 3 >      // Spline order
class ITK_EXPORT AdvancedBSplineDeformableTransform
  : public AdvancedTransform< TScalarType, NDimensions, NDimensions >
{
public:
  /** Standard class typedefs. */
  typedef AdvancedBSplineDeformableTransform        Self;
  typedef AdvancedTransform<
    TScalarType, NDimensions, NDimensions >         Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;
      
  /** New macro for creation of through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AdvancedBSplineDeformableTransform, AdvancedTransform );

  /** Dimension of the domain space. */
  itkStaticConstMacro( SpaceDimension, unsigned int, NDimensions );

  /** The BSpline order. */
  itkStaticConstMacro( SplineOrder, unsigned int, VSplineOrder );

  /** Typedefs from Superclass. */
  typedef typename Superclass::ParametersType         ParametersType;
  typedef typename Superclass::JacobianType           JacobianType;
  typedef typename Superclass::ScalarType             ScalarType;
  typedef typename Superclass::InputPointType         InputPointType;
  typedef typename Superclass::OutputPointType        OutputPointType;
  typedef typename Superclass::InputVectorType        InputVectorType;
  typedef typename Superclass::OutputVectorType       OutputVectorType;
  typedef typename Superclass::InputVnlVectorType     InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType    OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType 
    InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType
    OutputCovariantVectorType;

  typedef typename Superclass
    ::NonZeroJacobianIndicesType                    NonZeroJacobianIndicesType;
  typedef typename Superclass::SpatialJacobianType  SpatialJacobianType;
  typedef typename Superclass
    ::JacobianOfSpatialJacobianType                 JacobianOfSpatialJacobianType;
  typedef typename Superclass::SpatialHessianType   SpatialHessianType;
  typedef typename Superclass
    ::JacobianOfSpatialHessianType                  JacobianOfSpatialHessianType;
  typedef typename Superclass::InternalMatrixType   InternalMatrixType;
  
  /** This method sets the parameters of the transform.
   * For a BSpline deformation transform, the parameters are the BSpline 
   * coefficients on a sparse grid. 
   * 
   * The parameters are N number of N-D grid of coefficients. Each N-D grid 
   * is represented as a flat array of doubles 
   * (in the same configuration as an itk::Image).
   * The N arrays are then concatenated to form one parameter array.
   *
   * For efficiency, this transform does not make a copy of the parameters.
   * It only keeps a pointer to the input parameters. It assumes that the memory
   * is managed by the caller. Use SetParametersByValue to force the transform
   * to call copy the parameters.
   *
   * This method wraps each grid as itk::Image's using the user specified
   * grid region, spacing and origin.
   * NOTE: The grid region, spacing and origin must be set first.
   */
  void SetParameters( const ParametersType & parameters );
  
  /** This method sets the fixed parameters of the transform.
   * For a BSpline deformation transform, the parameters are the following:
   *    Grid Size, Grid Origin, and Grid Spacing
   * 
   * The fixed parameters are the three times the size of the templated 
   * dimensions.
   * This function has the effect of make the following calls:
   *       transform->SetGridSpacing( spacing );
   *       transform->SetGridOrigin( origin );
   *       transform->SetGridDirection( direction );
   *       transform->SetGridRegion( bsplineRegion );
   *
   * This function was added to allow the transform to work with the 
   * itkTransformReader/Writer I/O filters.
   */
  void SetFixedParameters( const ParametersType & parameters );

  /** This method sets the parameters of the transform.
   * For a BSpline deformation transform, the parameters are the BSpline 
   * coefficients on a sparse grid. 
   * 
   * The parameters are N number of N-D grid of coefficients. Each N-D grid 
   * is represented as a flat array of doubles 
   * (in the same configuration as an itk::Image).
   * The N arrays are then concatenated to form one parameter array.
   *
   * This methods makes a copy of the parameters while for
   * efficiency the SetParameters method does not.
   *
   * This method wraps each grid as itk::Image's using the user specified
   * grid region, spacing and origin.
   * NOTE: The grid region, spacing and origin must be set first.
   */
  void SetParametersByValue( const ParametersType & parameters );

  /** This method can ONLY be invoked AFTER calling SetParameters(). 
   *  This restriction is due to the fact that the AdvancedBSplineDeformableTransform
   *  does not copy the array of parameters internally, instead it keeps a 
   *  pointer to the user-provided array of parameters. This method is also
   *  in violation of the const-correctness of the parameters since the 
   *  parameter array has been passed to the transform on a 'const' basis but
   *  the values get modified when the user invokes SetIdentity().
   */
  void SetIdentity( void );

  /** Get the Transformation Parameters. */
  virtual const ParametersType& GetParameters( void ) const;

  /** Get the Transformation Fixed Parameters. */
  virtual const ParametersType& GetFixedParameters( void ) const;
  
  /** Parameters as SpaceDimension number of images. */
  typedef typename ParametersType::ValueType            PixelType;
  typedef Image< PixelType,
    itkGetStaticConstMacro( SpaceDimension )>           ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;

  /** Get the array of coefficient images. */
  virtual ImagePointer * GetCoefficientImage( void )
    { return this->m_CoefficientImage; }
  virtual const ImagePointer * GetCoefficientImage( void ) const
    { return this->m_CoefficientImage; }

  /** Set the array of coefficient images.
   *
   * This is an alternative API for setting the BSpline coefficients
   * as an array of SpaceDimension images. The grid region spacing 
   * and origin is taken from the first image. It is assume that
   * the buffered region of all the subsequent images are the same 
   * as the first image. Note that no error checking is done.
   *
   * Warning: use either the SetParameters() or SetCoefficientImage()
   * API. Mixing the two modes may results in unexpected results.
   */
  virtual void SetCoefficientImage( ImagePointer images[] );  

  /** Typedefs for specifying the extend to the grid. */
  typedef ImageRegion< itkGetStaticConstMacro( SpaceDimension ) > RegionType;
  
  typedef typename RegionType::IndexType    IndexType;
  typedef typename RegionType::SizeType     SizeType;
  typedef typename ImageType::SpacingType   SpacingType;
  typedef typename ImageType::DirectionType DirectionType;
  typedef typename ImageType::PointType     OriginType;
  typedef IndexType                         GridOffsetType;

  /** This method specifies the region over which the grid resides. */
  virtual void SetGridRegion( const RegionType& region );
  itkGetMacro( GridRegion, RegionType );
  itkGetConstMacro( GridRegion, RegionType );

  /** This method specifies the grid spacing or resolution. */
  virtual void SetGridSpacing( const SpacingType & spacing );
  itkGetMacro( GridSpacing, SpacingType );
  itkGetConstMacro( GridSpacing, SpacingType );

  /** This method specifies the grid directions . */
  virtual void SetGridDirection( const DirectionType & spacing );
  itkGetMacro( GridDirection, DirectionType );
  itkGetConstMacro( GridDirection, DirectionType );

  /** This method specifies the grid origin. */
  virtual void SetGridOrigin( const OriginType& origin );
  itkGetMacro( GridOrigin, OriginType );
  itkGetConstMacro( GridOrigin, OriginType );

  /** Typedef of the bulk transform. */
  typedef Transform< ScalarType,
    itkGetStaticConstMacro( SpaceDimension ),
    itkGetStaticConstMacro( SpaceDimension ) >          BulkTransformType;
  typedef typename BulkTransformType::ConstPointer      BulkTransformPointer;

  /** This method specifies the bulk transform to be applied. 
   * The default is the identity transform.
   */
  itkSetConstObjectMacro( BulkTransform, BulkTransformType );
  itkGetConstObjectMacro( BulkTransform, BulkTransformType );

  /** Transform points by a BSpline deformable transformation. */
  OutputPointType TransformPoint( const InputPointType & point ) const;

  /** Interpolation weights function type. */
  typedef BSplineInterpolationWeightFunction2< ScalarType,
    itkGetStaticConstMacro( SpaceDimension ),
    itkGetStaticConstMacro( SplineOrder ) >                 WeightsFunctionType;
  typedef typename WeightsFunctionType::WeightsType         WeightsType;
  typedef typename WeightsFunctionType::ContinuousIndexType ContinuousIndexType;
  typedef BSplineInterpolationDerivativeWeightFunction<
    ScalarType,
    itkGetStaticConstMacro( SpaceDimension ),
    itkGetStaticConstMacro( SplineOrder ) >                 DerivativeWeightsFunctionType;
  typedef BSplineInterpolationSecondOrderDerivativeWeightFunction<
    ScalarType,
    itkGetStaticConstMacro( SpaceDimension ),
    itkGetStaticConstMacro( SplineOrder ) >                 SODerivativeWeightsFunctionType;

  /** Parameter index array type. */
  typedef Array<unsigned long> ParameterIndexArrayType;

  /** Transform points by a BSpline deformable transformation. 
   * On return, weights contains the interpolation weights used to compute the 
   * deformation and indices of the x (zeroth) dimension coefficient parameters
   * in the support region used to compute the deformation.
   * Parameter indices for the i-th dimension can be obtained by adding
   * ( i * this->GetNumberOfParametersPerDimension() ) to the indices array.
   */
  virtual void TransformPoint(
    const InputPointType & inputPoint,
    OutputPointType & outputPoint,
    WeightsType & weights,
    ParameterIndexArrayType & indices,
    bool & inside ) const;

  /** Get number of weights. */
  unsigned long GetNumberOfWeights( void ) const
  {
    return this->m_WeightsFunction->GetNumberOfWeights();    
  }

  /** Method to transform a vector - 
   *  not applicable for this type of transform.
   */
  virtual OutputVectorType TransformVector( const InputVectorType & ) const
    { 
    itkExceptionMacro( << "Method not applicable for deformable transform." );
    return OutputVectorType(); 
    }

  /** Method to transform a vnl_vector - 
   *  not applicable for this type of transform.
   */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType & ) const
    { 
    itkExceptionMacro( << "Method not applicable for deformable transform. ");
    return OutputVnlVectorType(); 
    }

  /** Method to transform a CovariantVector - 
   *  not applicable for this type of transform.
   */
  virtual OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType & ) const
    { 
    itkExceptionMacro( << "Method not applicable for deformable transform. ");
    return OutputCovariantVectorType(); 
    } 
    
  /** Return the number of parameters that completely define the Transform. */
  virtual unsigned int GetNumberOfParameters( void ) const;

  /** Return the number of parameters per dimension */
  unsigned int GetNumberOfParametersPerDimension( void ) const;

  /** Return the region of the grid wholly within the support region */
  itkGetConstReferenceMacro( ValidRegion, RegionType );

  /** Indicates that this transform is linear. That is, given two
   * points P and Q, and scalar coefficients a and b, then
   *
   *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
   */
  virtual bool IsLinear( void ) const { return false; }

  unsigned int GetNumberOfAffectedWeights( void ) const;

  virtual unsigned long GetNumberOfNonZeroJacobianIndices( void ) const;

  /** Whether the advanced transform has nonzero matrices. */
  virtual bool GetHasNonZeroSpatialHessian( void ) const
  {
    return true;
  }
  virtual bool HasNonZeroJacobianOfSpatialHessian( void ) const
  {
    return true;
  }

  /** Compute the Jacobian matrix of the transformation at one point. */
  virtual const JacobianType & GetJacobian( const InputPointType & point ) const;

  /** Compute the Jacobian of the transformation. */
  virtual void GetJacobian(
    const InputPointType & ipp,
    WeightsType & weights,
    ParameterIndexArrayType & indices ) const;

  /** Compute the Jacobian of the transformation. */
  virtual void GetJacobian(
    const InputPointType & ipp,
    JacobianType & j,
    NonZeroJacobianIndicesType & ) const;

  /** Compute the spatial Jacobian of the transformation. */
  virtual void GetSpatialJacobian(
    const InputPointType & ipp,
    SpatialJacobianType & sj ) const;

  /** Compute the spatial Hessian of the transformation. */
  virtual void GetSpatialHessian(
    const InputPointType & ipp,
    SpatialHessianType & sh ) const;

  /** Compute the Jacobian of the spatial Jacobian of the transformation. */
  virtual void GetJacobianOfSpatialJacobian(
    const InputPointType & ipp,
    JacobianOfSpatialJacobianType & jsj,
    NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const;

  /** Compute both the spatial Jacobian and the Jacobian of the
   * spatial Jacobian of the transformation.
   */
  virtual void GetJacobianOfSpatialJacobian(
    const InputPointType & ipp,
    SpatialJacobianType & sj,
    JacobianOfSpatialJacobianType & jsj,
    NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const;

  /** Compute the Jacobian of the spatial Hessian of the transformation. */
  virtual void GetJacobianOfSpatialHessian(
    const InputPointType & ipp,
    JacobianOfSpatialHessianType & jsh,
    NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const;

  /** Compute both the spatial Hessian and the Jacobian of the
   * spatial Hessian of the transformation.
   */
  virtual void GetJacobianOfSpatialHessian(
    const InputPointType & ipp,
    SpatialHessianType & sh,
    JacobianOfSpatialHessianType & jsh,
    NonZeroJacobianIndicesType & nonZeroJacobianIndices ) const;

protected:
  /** Print contents of an AdvancedBSplineDeformableTransform. */
  virtual void PrintSelf( std::ostream &os, Indent indent ) const;

  AdvancedBSplineDeformableTransform();
  virtual ~AdvancedBSplineDeformableTransform();

  /** Allow subclasses to access and manipulate the weights function. */
  // Why??
  itkSetObjectMacro( WeightsFunction, WeightsFunctionType );
  itkGetObjectMacro( WeightsFunction, WeightsFunctionType );

  /** Wrap flat array into images of coefficients. */
  void WrapAsImages( void );

  /** Convert an input point to a continuous index inside the BSpline grid. */
  void TransformPointToContinuousGridIndex( 
   const InputPointType & point, ContinuousIndexType & index ) const;

  virtual void ComputeNonZeroJacobianIndices(
    NonZeroJacobianIndicesType & nonZeroJacobianIndices,
    const RegionType & supportRegion ) const;

  /** Check if a continuous index is inside the valid region. */
  virtual bool InsideValidRegion( const ContinuousIndexType& index ) const;

  /** The bulk transform. */
  BulkTransformPointer  m_BulkTransform;

  /** Array of images representing the B-spline coefficients 
  *  in each dimension. */
  ImagePointer    m_CoefficientImage[ NDimensions ];

  /** Pointer to function used to compute B-spline interpolation weights. */
  typename WeightsFunctionType::Pointer             m_WeightsFunction;
  typename DerivativeWeightsFunctionType::Pointer   m_DerivativeWeightsFunction;
  typename SODerivativeWeightsFunctionType::Pointer m_SODerivativeWeightsFunction;

  /** Variables defining the coefficient grid extend. */
  RegionType      m_GridRegion;
  SpacingType     m_GridSpacing;
  DirectionType   m_GridDirection;
  OriginType      m_GridOrigin;
  GridOffsetType  m_GridOffsetTable;

  DirectionType   m_PointToIndexMatrix;
  DirectionType   m_IndexToPoint;

  RegionType      m_ValidRegion;

  /** Variables defining the interpolation support region. */
  unsigned long   m_Offset;
  bool            m_SplineOrderOdd;
  SizeType        m_SupportSize;
  IndexType       m_ValidRegionLast;

  /** Keep a pointer to the input parameters. */
  const ParametersType *  m_InputParametersPointer;

  /** Jacobian as SpaceDimension number of images. */
  typedef typename JacobianType::ValueType      JacobianPixelType;
  typedef Image< JacobianPixelType,
    itkGetStaticConstMacro( SpaceDimension ) >  JacobianImageType;

  typename JacobianImageType::Pointer m_JacobianImage[ NDimensions ];

  /** Keep track of last support region used in computing the Jacobian
  * for fast resetting of Jacobian to zero.
  */
  mutable IndexType m_LastJacobianIndex;

  /** Array holding images wrapped from the flat parameters. */
  ImagePointer    m_WrappedImage[ NDimensions ];

  /** Internal parameters buffer. */
  ParametersType          m_InternalParametersBuffer;

private:
  AdvancedBSplineDeformableTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  void UpdateGridOffsetTable( void );

}; //class AdvancedBSplineDeformableTransform


}  // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_AdvancedBSplineDeformableTransform(_, EXPORT, x, y) namespace itk { \
  _(3(class EXPORT AdvancedBSplineDeformableTransform< ITK_TEMPLATE_3 x >)) \
  namespace Templates { typedef AdvancedBSplineDeformableTransform< ITK_TEMPLATE_3 x > \
                                                   AdvancedBSplineDeformableTransform##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAdvancedBSplineDeformableTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAdvancedBSplineDeformableTransform.txx"
#endif

#endif /* __itkAdvancedBSplineDeformableTransform_h */
