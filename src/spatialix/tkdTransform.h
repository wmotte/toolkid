#ifndef __tkdTransform_h__
#define __tkdTransform_h__
#include "itkAdvancedTransform.h"
#include "itkParameterMapInterface.h"
//#include "itkImage.h"

namespace tkd
{

class Transform
  : public itk::LightObject
{
public:
  typedef Transform Self;
  typedef itk::LightObject Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkTypeMacro( Transform, LightObject );

  static const unsigned int Dimension = 3;
  typedef double ScalarType;

  typedef itk::AdvancedTransform< ScalarType, Dimension, Dimension > TransformType;
  typedef TransformType::Pointer TransformPointer;
  typedef TransformType::ConstPointer TransformConstPointer;
  typedef TransformType::ParametersType TransformParametersType;

  typedef std::vector< std::string > ParameterValuesType;
  typedef std::map< std::string, ParameterValuesType > ParameterMapType;

  virtual void SetParameterMap( const ParameterMapType& map );
  virtual const ParameterMapType& GetParameterMap();

  TransformPointer GetTransform();
  virtual Pointer GetInverse();
  virtual void Initialize( itk::ParameterMapInterface::Pointer parameters );
  virtual void Log( std::ostream& os );
  virtual int Write( const std::string& path, int counter = 0, bool compose = true );
  virtual void Write( std::ostream& os );

protected:
  Transform();
  virtual ~Transform();

  virtual void Clone( Self* pointer );

  TransformPointer m_Transform;
  ParameterMapType m_ParameterMap;
  TransformParametersType m_TransformParameters;

  template< class T >
  void WriteParameter( const std::string& parameterName, const T& value, std::ostream& os, bool string = true )
  {
    os << "(" << parameterName << " " << ( string ? "\"" : "" ) << value << ( string ? "\"" : "" ) << ")" << std::endl;
  }

  template< class T >
  void WriteParameters( const std::string& parameterName, T value, int length, std::ostream& os, bool string = true )
  {
    os << "(" << parameterName;
    for( int i = 0; i < length; ++i )
      {
      os << " ";
      os << ( string ? "\"" : "" );
      os << value[ i ];
      os << ( string ? "\"" : "" );
      }
    os << ")" << std::endl;
  }

private:
  Transform( const Self& );
  void operator=( const Self& );
};

} // end namespace tkd

#endif /*__tkdTransform_h__*/
