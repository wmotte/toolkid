#ifndef __tkdBSplineTransform_h__
#define __tkdBSplineTransform_h__
#include "tkdTransform.h"
#include "itkAdvancedBSplineDeformableTransform.h"

namespace tkd
{

class BSplineTransform : public Transform
{
public:
  typedef BSplineTransform Self;
  typedef Transform Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( BSplineTransform, Transform );

  typedef Superclass::TransformType TransformType;
  typedef Superclass::ScalarType ScalarType;
  static const unsigned int Dimension = Superclass::Dimension;

  typedef itk::AdvancedBSplineDeformableTransform< ScalarType, Dimension, Dimension > BSplineTransformType;
  typedef BSplineTransformType::Pointer BSplineTransformPointer;
  typedef BSplineTransformType::ConstPointer BSplineTransformConstPointer;

  virtual void Initialize( itk::ParameterMapInterface::Pointer parameters );
  virtual void Log( std::ostream& os );
  virtual void Write( std::ostream& os );

protected:
  BSplineTransform();
  virtual ~BSplineTransform();

  BSplineTransformPointer m_BSplineTransform;
  BSplineTransformType::ParametersType m_Parameters;

private:
  BSplineTransform( const Self& );
  void operator=( const Self& );
};

} // end namespace tkd

#endif /*__tkdBSplineTransform_h__*/
