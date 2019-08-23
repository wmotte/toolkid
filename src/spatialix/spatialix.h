#ifndef __tkdSpatialix_h__
#define __tkdSpatialix_h__
#include "tkdTransform.h"
#include "tkdTransformFactory.h"

namespace tkd
{

class Spatialix
{
public:
  Spatialix();

  Transform::Pointer Read( const std::vector< std::string >& transformFileNames );

  /**
   * Read transform from file
   * @param readInitial Also load initial transform (true)
   * @param allowMatrixComposition Allow matrix composition of two matrix transforms are combined by composition
   */
  Transform::Pointer Read(
      const std::string& transformFileName,
      bool readInitial = true,
      bool allowMatrixComposition = true );

  /**
   * Combine two transforms.
   * @param compose Combination by composition (true) or addition (false)
   * @param allowMatrixComposition Allow matrix composition of two matrix transforms are combined by composition
   */
  Transform::Pointer Combine(
      Transform::Pointer initial,
      Transform::Pointer current,
      bool compose = true,
      bool allowMatrixComposition = true );

  Transform::Pointer NormalizeBSplineTransform( Transform::Pointer transform );

protected:
  TransformFactory::Pointer m_Factory;
};

} // end namespace tkd

#endif /*__tkdSpatialix_h__*/
