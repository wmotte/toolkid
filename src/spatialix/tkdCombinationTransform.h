#ifndef __tkdCombinationTransform_h__
#define __tkdCombinationTransform_h__
#include "tkdTransform.h"
#include "itkAdvancedCombinationTransform.h"

namespace tkd
{

class CombinationTransform : public Transform
{
public:
  typedef CombinationTransform Self;
  typedef Transform Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( CombinationTransform, Transform );

  typedef Superclass::TransformType TransformType;
  typedef Superclass::ScalarType ScalarType;
  static const unsigned int Dimension = Superclass::Dimension;

  typedef itk::AdvancedCombinationTransform< ScalarType, Dimension > CombinationTransformType;
  typedef CombinationTransformType::Pointer CombinationTransformPointer;
  typedef CombinationTransformType::ConstPointer CombinationTransformConstPointer;

  void SetInitialTransform( Superclass::Pointer transform );
  void SetCurrentTransform( Superclass::Pointer transform );

  Superclass::Pointer GetInitialTransform();
  Superclass::Pointer GetCurrentTransform();

  bool GetUseComposition();
  void SetUseComposition( bool composition );

  CombinationTransformPointer GetCombinationTransform();

  virtual Superclass::Pointer GetInverse();
  virtual void Initialize( itk::ParameterMapInterface::Pointer parameters );
  virtual void Log( std::ostream& os );
  virtual int Write( const std::string& path, int counter = 0, bool compose = true );

protected:
  CombinationTransform();
  virtual ~CombinationTransform();

  CombinationTransformPointer m_CombinationTransform;
  Superclass::Pointer m_InitialTransform;
  Superclass::Pointer m_CurrentTransform;

private:
  CombinationTransform( const Self& );
  void operator=( const Self& );
};

} // end namespace tkd

#endif /*__tkdCombinationTransform_h__*/
