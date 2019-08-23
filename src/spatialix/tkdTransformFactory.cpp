#include "tkdTransformFactory.h"
#include "tkdBSplineTransform.h"
#include "tkdSimilarityTransform.h"
#include "tkdAffineTransform.h"

namespace tkd
{

TransformFactory::TransformFactory()
{
}

TransformFactory::~TransformFactory()
{
}

Transform::Pointer TransformFactory::CreateTransform( const std::string& transformName )
{
  if ( transformName == "BSplineTransform" )
    {
    return tkd::BSplineTransform::New().GetPointer();
    }

  if ( transformName == "SimilarityTransform" )
    {
    return tkd::SimilarityTransform::New().GetPointer();
    }

  if ( transformName == "AffineTransform" )
    {
    return tkd::AffineTransform::New().GetPointer();
    }

  return 0;
}

} // end namespace tkd
