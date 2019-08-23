#include "tkdBSplineTransform.h"

namespace tkd
{

BSplineTransform::BSplineTransform()
{
  this->m_BSplineTransform = BSplineTransformType::New();
  this->m_Transform = m_BSplineTransform.GetPointer();
}

BSplineTransform::~BSplineTransform()
{

}

void BSplineTransform::Initialize( itk::ParameterMapInterface::Pointer parameters )
{
  BSplineTransformType::RegionType region;
  BSplineTransformType::IndexType index;
  BSplineTransformType::SizeType size;
  BSplineTransformType::SpacingType spacing;
  BSplineTransformType::OriginType origin;

  std::string error;

  for( unsigned int i = 0; i < Dimension; ++i )
    {
    parameters->ReadParameter< double >( origin[ i ], "GridOrigin", i, true, error );
    parameters->ReadParameter< unsigned long >( size[ i ], "GridSize", i, true, error );
    parameters->ReadParameter< double >( spacing[ i ], "GridSpacing", i, true, error );
    parameters->ReadParameter< long >( index[ i ], "GridIndex", i, true, error );
    }

  region.SetIndex( index );
  region.SetSize( size );

  this->m_BSplineTransform->SetGridRegion( region );
  this->m_BSplineTransform->SetGridOrigin( origin );
  this->m_BSplineTransform->SetGridSpacing( spacing );

  Superclass::Initialize( parameters );
}

void BSplineTransform::Log( std::ostream& os )
{
  Superclass::Log( os );
  this->m_BSplineTransform->Print( os );
}

void BSplineTransform::Write( std::ostream& os )
{
  this->WriteParameter< std::string >( "Transform", "BSplineTransform", os, true );
  this->WriteParameters< BSplineTransformType::SizeType >( "GridSize", this->m_BSplineTransform->GetGridRegion().GetSize(), Dimension, os, false );
  this->WriteParameters< BSplineTransformType::IndexType >( "GridIndex", this->m_BSplineTransform->GetGridRegion().GetIndex(), Dimension, os, false );
  this->WriteParameters< BSplineTransformType::OriginType >( "GridOrigin", this->m_BSplineTransform->GetGridOrigin(), Dimension, os, false );
  this->WriteParameters< BSplineTransformType::SpacingType >( "GridSpacing", this->m_BSplineTransform->GetGridSpacing(), Dimension, os, false );

  Superclass::Write( os );
}

} // end namespace tkd
