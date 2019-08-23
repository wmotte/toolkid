#include "tkdSimilarityTransform.h"

namespace tkd
{

SimilarityTransform::SimilarityTransform()
{
  this->m_SimilarityTransform = SimilarityTransformType::New();
  this->m_Transform = this->m_SimilarityTransform.GetPointer();
}

SimilarityTransform::~SimilarityTransform()
{

}

SimilarityTransform::Superclass::Pointer SimilarityTransform::GetInverse()
{
  SimilarityTransformType::Pointer inverse = SimilarityTransformType::New();
  inverse->SetCenter( this->m_SimilarityTransform->GetCenter() );
  if ( this->m_SimilarityTransform->GetInverse( inverse ) )
    {
    Pointer transform = New();
    transform->Clone( this );
    transform->m_SimilarityTransform = inverse;
    transform->m_Transform = inverse.GetPointer();
    return transform.GetPointer();
    }

  return 0;
}

void SimilarityTransform::Initialize( itk::ParameterMapInterface::Pointer parameters )
{
  SimilarityTransformType::InputPointType center;
  std::string error;

  for( unsigned int i = 0; i < Dimension; ++i )
    {
    parameters->ReadParameter< double >( center[ i ], "CenterOfRotationPoint", i, true, error );
    }

  this->m_SimilarityTransform->SetCenter( center );

  Superclass::Initialize( parameters );
}

void SimilarityTransform::Log( std::ostream& os )
{
  Superclass::Log( os );
  this->m_SimilarityTransform->Print( os );
}

void SimilarityTransform::Write( std::ostream& os )
{
  this->WriteParameter< std::string >( "Transform", "SimilarityTransform", os, true );
  this->WriteParameters< SimilarityTransformType::InputPointType >( "CenterOfRotationPoint", this->m_SimilarityTransform->GetCenter(), Dimension, os, false );
  Superclass::Write( os );
}

} // end namespace tkd
