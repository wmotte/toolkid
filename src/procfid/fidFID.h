#ifndef __fidFID_h__
#define __fidFID_h__
#include "itkObjectFactory.h"
#include "itkLightObject.h"
#include "itkAutoPointer.h"
#include "itkImage.h"
#include <complex>
#include "vnl/vnl_matrix_ref.h"
#include "procparser.h"

namespace fid
{

class FIDImpl;

class FID : public itk::LightObject
{
public:
  typedef FID Self;
  typedef itk::LightObject Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( FID, LightObject );

  typedef float PrecisionType;
  typedef std::complex< PrecisionType > ComplexType;
  typedef itk::Image< ComplexType, 3 > DataType;
  typedef DataType::Pointer DataPointer;
  typedef DataType::ConstPointer DataConstPointer;
  typedef vnl_matrix_ref< ComplexType > BlockType;

  BlockType GetBlock( int block );

  DataPointer GetData();
  DataConstPointer GetData() const;

  int GetNumberOfBlocks() const;
  int GetNumberOfTraces() const;
  int GetNumberOfPoints() const;

  const Procparser& GetProcpar() const;

  template< class T >
  T Procpar( const std::string& key, int i = 0, T d = T() ) const
  {
    const Procparser& pp = GetProcpar();
    return pp.GetAs< T >( key, i, d );
  }

protected:
  FID();
  virtual ~FID();

  DataPointer m_Data;

private:
  FID( const Self& );
  void operator=( const Self& );

  friend class FIDReader;
  itk::AutoPointer< FIDImpl > m_Impl;
};

} // end namespace fid

#endif /*__fidFID_h__*/
