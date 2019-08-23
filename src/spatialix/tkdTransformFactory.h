#ifndef __tkdTransformFactory_h__
#define __tkdTransformFactory_h__
#include "tkdTransform.h"

namespace tkd
{

class TransformFactory
  : public itk::LightObject
{
public:
  typedef TransformFactory Self;
  typedef itk::LightObject Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( TransformFactory, LightObject );

  Transform::Pointer CreateTransform( const std::string& transformName );

protected:
  TransformFactory();
  virtual ~TransformFactory();

private:
  TransformFactory( const Self& );
  void operator=( const Self& );
};

} // end namespace tkd

#endif /*__tkdTransformFactory_h__*/
