#include "tkdAffineTransform.h"

namespace tkd
{

AffineTransform::AffineTransform()
{
  this->m_AffineTransform = AffineTransformType::New();
  this->m_Transform = this->m_AffineTransform.GetPointer();
  this->m_AffineTransform->SetIdentity();
}

AffineTransform::~AffineTransform()
{

}

AffineTransform::Superclass::Pointer AffineTransform::GetInverse()
{
  AffineTransformType::Pointer inverse = AffineTransformType::New();
  inverse->SetCenter( this->m_AffineTransform->GetCenter() );
  if ( this->m_AffineTransform->GetInverse( inverse ) )
    {
    Pointer transform = New();
    transform->Clone( this );
    transform->m_AffineTransform = inverse;
    transform->m_Transform = inverse.GetPointer();
    return transform.GetPointer();
    }

  return 0;
}

void AffineTransform::Initialize( itk::ParameterMapInterface::Pointer parameters )
{
  AffineTransformType::InputPointType center;
  std::string error;

  for( unsigned int i = 0; i < Dimension; ++i )
    {
    parameters->ReadParameter< double >( center[ i ], "CenterOfRotationPoint", i, true, error );
    }

  this->m_AffineTransform->SetCenter( center );

  Superclass::Initialize( parameters );
}

void AffineTransform::Log( std::ostream& os )
{
  Superclass::Log( os );
  this->m_AffineTransform->Print( os );
}

void AffineTransform::Write( std::ostream& os )
{
  this->WriteParameter< std::string >( "Transform", "AffineTransform", os, true );
  this->WriteParameters< AffineTransformType::InputPointType >( "CenterOfRotationPoint", this->m_AffineTransform->GetCenter(), Dimension, os, false );
  Superclass::Write( os );
}

AffineTransform::AffineTransformPointer AffineTransform::GetAffineTransform()
{
  return this->m_AffineTransform;
}

} // end namespace tkd
