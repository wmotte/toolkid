#ifndef __tkdSimilarityTransform_h__
#define __tkdSimilarityTransform_h__
#include "tkdTransform.h"
#include "itkAdvancedSimilarity3DTransform.h"

namespace tkd
{

class SimilarityTransform : public Transform
{
public:
  typedef SimilarityTransform Self;
  typedef Transform Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( SimilarityTransform, Transform );

  typedef Superclass::TransformType TransformType;
  typedef Superclass::ScalarType ScalarType;
  static const unsigned int Dimension = Superclass::Dimension;

  typedef itk::AdvancedSimilarity3DTransform< ScalarType > SimilarityTransformType;
  typedef SimilarityTransformType::Pointer SimilarityTransformPointer;
  typedef SimilarityTransformType::ConstPointer SimilarityTransformConstPointer;

  virtual Superclass::Pointer GetInverse();
  virtual void Initialize( itk::ParameterMapInterface::Pointer parameters );
  virtual void Log( std::ostream& os );
  virtual void Write( std::ostream& os );

protected:
  SimilarityTransform();
  virtual ~SimilarityTransform();

  SimilarityTransformPointer m_SimilarityTransform;

private:
  SimilarityTransform( const Self& );
  void operator=( const Self& );
};

} // end namespace tkd

#endif /*__tkdSimilarityTransform_h__*/
