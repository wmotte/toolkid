#include "tkdCombinationTransform.h"

namespace tkd
{

CombinationTransform::CombinationTransform()
{
  this->m_CombinationTransform = CombinationTransformType::New();
  this->m_Transform = this->m_CombinationTransform.GetPointer();
}

CombinationTransform::~CombinationTransform()
{

}

void CombinationTransform::SetUseComposition( bool composition )
{
  this->m_CombinationTransform->SetUseComposition( composition );
}

bool CombinationTransform::GetUseComposition()
{
  return this->m_CombinationTransform->GetUseComposition();
}

CombinationTransform::CombinationTransformPointer CombinationTransform::GetCombinationTransform()
{
  return this->m_CombinationTransform;
}

CombinationTransform::Superclass::Pointer CombinationTransform::GetInitialTransform()
{
  return this->m_InitialTransform;
}

CombinationTransform::Superclass::Pointer CombinationTransform::GetCurrentTransform()
{
  return this->m_CurrentTransform;
}

CombinationTransform::Superclass::Pointer CombinationTransform::GetInverse()
{
  if ( !this->m_CurrentTransform )
    {
    return 0;
    }

  if ( !this->m_InitialTransform )
    {
    return this->m_CurrentTransform->GetInverse();
    }

  if( !this->GetUseComposition() )
    {
    return 0;
    }

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
  Superclass::Pointer inverseT0 = this->m_InitialTransform->GetInverse();

  if ( inverseT0 )
    {
    /** Try to create the inverse of the current transform. */
    Superclass::Pointer inverseT1 = this->m_CurrentTransform->GetInverse();

    if ( inverseT1 )
      {
      Pointer transform = New();
      transform->SetInitialTransform( inverseT1 );
      transform->SetCurrentTransform( inverseT0 );
      transform->SetUseComposition( true );
      return transform.GetPointer();
      }
    }

  return 0;
}

void CombinationTransform::Initialize( itk::ParameterMapInterface::Pointer parameters )
{

}

void CombinationTransform::SetInitialTransform( Superclass::Pointer transform )
{
  this->m_InitialTransform = transform;
  this->m_CombinationTransform->SetInitialTransform( transform->GetTransform() );
}

void CombinationTransform::SetCurrentTransform( Superclass::Pointer transform )
{
  this->m_CurrentTransform = transform;
  this->m_CombinationTransform->SetCurrentTransform( transform->GetTransform() );
}

void CombinationTransform::Log( std::ostream& os )
{
  Superclass::Log( os );
  os << "[" << this << "] Combination method: " << ( this->GetUseComposition() ? "Compose" : "Add" ) << std::endl;
  os << "[" << this << "] Current transform:" << std::endl;
  this->m_CurrentTransform->Log( os );
  os << std::endl;
  os << "[" << this << "] Initial transform:" << std::endl;
  this->m_InitialTransform->Log( os );
  os << std::endl;
}

int CombinationTransform::Write( const std::string& path, int counter, bool compose )
{
  if ( this->m_InitialTransform )
    {
    counter = this->m_InitialTransform->Write( path, counter );
    }

  return this->m_CurrentTransform->Write( path, counter, this->GetUseComposition() );
}

} // end namespace tkd
