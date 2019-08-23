#ifndef __tkdAffineTransform_h__
#define __tkdAffineTransform_h__
#include "tkdTransform.h"
#include "itkAdvancedMatrixOffsetTransformBase.h"

namespace tkd
{

class AffineTransform : public Transform
{
public:
  typedef AffineTransform Self;
  typedef Transform Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( AffineTransform, Transform );

  typedef Superclass::TransformType TransformType;
  typedef Superclass::ScalarType ScalarType;
  static const unsigned int Dimension = Superclass::Dimension;

  typedef itk::AdvancedMatrixOffsetTransformBase< ScalarType, Dimension, Dimension > AffineTransformType;
  typedef AffineTransformType::Pointer AffineTransformPointer;
  typedef AffineTransformType::ConstPointer AffineTransformConstPointer;

  virtual void Initialize( itk::ParameterMapInterface::Pointer parameters );

  virtual Superclass::Pointer GetInverse();
  virtual void Log( std::ostream& os );
  virtual void Write( std::ostream& os );

  AffineTransformPointer GetAffineTransform();

protected:
  AffineTransform();
  virtual ~AffineTransform();

  AffineTransformPointer m_AffineTransform;

private:
  AffineTransform( const Self& );
  void operator=( const Self& );
};

} // end namespace tkd

#endif /*__tkdAffineTransform_h__*/
