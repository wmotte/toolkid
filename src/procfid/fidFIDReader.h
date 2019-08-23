#ifndef __fidFIDReader_h__
#define __fidFIDReader_h__
#include "itkObjectFactory.h"
#include "fidFID.h"

namespace fid
{

class FIDReader : public itk::LightObject
{
public:
  typedef FIDReader Self;
  typedef itk::LightObject Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( FIDReader, LightObject );

  typedef FID::PrecisionType PrecisionType;
  typedef FID::ComplexType ComplexType;

  void SetFileName( const std::string& filename );
  const std::string& GetFileName() const;

  void Read();
  FID::Pointer GetFID();

protected:
  FIDReader();
  virtual ~FIDReader();

  template< class TPrecisionInput >
  void ReadFID( std::ifstream& in, fid::FID::Pointer fid );

  std::string GetFullPath( const std::string& file );

  std::string m_FileName;
  FID::Pointer m_FID;

private:
  FIDReader( const Self& );
  void operator=( const Self& );
};

} // end namespace fid

#endif /*__fidFIDReader_h__*/
