#include "tkdTransform.h"

namespace tkd
{

Transform::Transform()
{
}

Transform::~Transform()
{

}

Transform::TransformPointer Transform::GetTransform()
{
  return this->m_Transform;
}

Transform::Pointer Transform::GetInverse()
{
  return 0;
}

void Transform::SetParameterMap( const ParameterMapType& map )
{
  this->m_ParameterMap = map;
}

const Transform::ParameterMapType& Transform::GetParameterMap()
{
  return this->m_ParameterMap;
}

void Transform::Initialize( itk::ParameterMapInterface::Pointer parameters )
{
  std::string error;

  int numberOfParameters = 0;
  if ( !parameters->ReadParameter< int >(
          numberOfParameters,
          "NumberOfParameters",
          0, true, error ) )
    {
    std::cout << error << std::endl;
    }

  this->m_TransformParameters = TransformParametersType( numberOfParameters );
  for( int i = 0; i < numberOfParameters; ++i )
    {
    parameters->ReadParameter< TransformType::ScalarType >(
        this->m_TransformParameters[ i ],
        "TransformParameters",
        i,
        true, error );
    }

  this->m_Transform->SetParameters( this->m_TransformParameters );
}

void Transform::Log( std::ostream& os )
{
  os << "[" << this << "] " << this->GetNameOfClass() << std::endl;
}

int Transform::Write( const std::string& path, int counter, bool compose )
{
  std::stringstream ss;
  std::stringstream ssInitial;

  ss << path << "/TransformParameters." << counter << ".txt";

  std::ofstream out( ss.str().c_str() );

  this->Write( out );

  if ( counter > 0 )
    {
    ssInitial << path << "/TransformParameters." << ( counter - 1 ) << ".txt";
    }
  else
    {
    ssInitial << "NoInitialTransform";
    }

  this->WriteParameter< std::string >( "InitialTransformParametersFileName", ssInitial.str(), out, true );
  this->WriteParameter< std::string >( "HowToCombineTransforms", ( compose ? "Compose" : "Add" ), out, true );

  return counter + 1;
}

void Transform::Write( std::ostream& os )
{
  this->WriteParameter< int >( "NumberOfParameters", this->m_Transform->GetParameters().Size(), os, false );
  this->WriteParameters< TransformType::ParametersType >( "TransformParameters", this->m_Transform->GetParameters(), this->m_Transform->GetParameters().Size(), os, false );
}

void Transform::Clone( Self* pointer )
{

}

} // end namespace tkd
