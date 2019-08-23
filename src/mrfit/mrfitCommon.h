#ifndef __mrfitCommon_h__
#define __mrfitCommon_h__
#include "itkImage.h"
#include "itkLightObject.h"
#include "itkObjectFactory.h"
#include "itkVectorImage.h"
#include <vector>
#include "itkAutoPointer.h"

namespace mrfit
{

class MRFit : public itk::LightObject
{
public:
  typedef float PixelType;
  static const unsigned int InputDimension = 4;
  static const unsigned int MapDimension = 3;
  typedef itk::Image< PixelType, InputDimension > InputImageType;
  typedef itk::Image< PixelType, MapDimension > MapImageType;
  typedef InputImageType::Pointer InputImagePointer;
  typedef MapImageType::Pointer MapImagePointer;
  typedef std::vector< double > TimesArray;

  typedef MRFit Self;
  typedef itk::LightObject Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( MRFit, LightObject );

  void SetInput( InputImagePointer input );
  void SetTimes( const TimesArray& times );

  enum T1FittingType
  {
    IDEAL_STEADY_STATE = 0,
    INVERSION_RECOVERY = 1,
    ABSOLUTE_INVERSION_RECOVERY = 2,
    LOOK_LOCKER = 3,
    ABSOLUTE_LOOK_LOCKER = 4,
    HYBRID_STEADY_STATE_3PARAM = 5,
    INVERSION_RECOVERY_3PARAM = 6,
    ABSOLUTE_INVERSION_RECOVERY_3PARAM = 7
  };

  enum T2FittingType
  {
    LINEAR = 0,
    NON_LINEAR = 1,
    NON_LINEAR_WITH_CONSTANT = 2
  };

  enum T1OutputType
  {
	  T1_MAP 		= 0,
	  A_MAP 		= 1,
	  B_MAP 		= 2,
	  R_SQUARED_MAP = 3

  };

  void FitT1( T1FittingType type = LOOK_LOCKER, double maxT1 = 10.0 );
  void FitT2( T2FittingType type = LINEAR, double maxT2 = 10.0 );

  MapImagePointer GetMap( int index );
  int GetNumberOfMaps() const;

protected:
  MRFit();
  virtual ~MRFit();

private:
  struct Impl;
  itk::AutoPointer< Impl > m_Impl;

  MRFit( const Self& );
  void operator=( const Self& );
};

} // end namespace mrfit


#endif /*__mrfitCommon_h__*/
