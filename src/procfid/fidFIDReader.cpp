#include "fidFIDReader.h"
#include "fidFIDImpl.h"
#include <fstream>
#include "vnl/vnl_matrix.h"
#include <complex>
#include "itksys/SystemTools.hxx"
#include "fidCommon.h"

inline void ByteSwap( char* block, int size )
{
  register int i = 0;
  register int j = size - 1;

  for( ; i < j ; ++i, --j )
  {
    std::swap( block[ i ], block[ j ] );
  }
}

namespace fid
{

FIDReader::FIDReader()
{

}

FIDReader::~FIDReader()
{

}

std::string FIDReader::GetFullPath( const std::string& file )
{
  std::string base = this->GetFileName();
  if ( base.length() > 1 && base.substr( base.length() - 1, 1 ) == "/" )
    {
    base = base.substr( 0, base.length() - 1 );
    }

  if ( base.length() > 3 && base.substr( base.length() - 4, 4 ) == ".fid" )
    {
    base = base.substr( 0, base.length() - 4 );
    }

  if ( base.length() > 2 && base.substr( base.length() - 3, 3 ) == "fid"
  		 && itksys::SystemTools::FileExists( base.c_str(), true ) )
  	{
  	base = base.substr( 0, base.length() - 3 );
  	}
  else
  	{
  	if ( !itksys::SystemTools::FileExists( base.c_str(), false ) &&
  			 !itksys::SystemTools::FileExists( base.c_str(), true ) )
  		{
  		std::stringstream ss;
  		ss << base << ".fid";
  		base = ss.str();
  		}

  	}

  std::stringstream ss;
  ss << base << "/" << file;

  return ss.str();
}

void FIDReader::SetFileName( const std::string& filename )
{
  m_FileName = filename;
}

const std::string& FIDReader::GetFileName() const
{
  return m_FileName;
}

template< class TPrecisionInput >
void
FIDReader::ReadFID( std::ifstream& in, fid::FID::Pointer fid )
{
  // prepare data matrix
  typedef TPrecisionInput InputPrecisionType;
  typedef typename fid::FID::PrecisionType OutputPrecisionType;
  typedef typename fid::FID::ComplexType ComplexType;
  typedef typename fid::FID::DataType DataType;

  datafilehead& header = fid->m_Impl->Header;

  InputPrecisionType* buffer = new InputPrecisionType[ header.np ];

  fid->m_Data = DataType::New();
  typename DataType::RegionType region;
  typename DataType::SizeType size;

  size[ 0 ] = header.np / 2;
  size[ 1 ] = header.ntraces;
  size[ 2 ] = header.nblocks;
  region.SetSize( size );

  fid->m_Data->SetRegions( region );
  fid->m_Data->Allocate();

  ComplexType* dataPointer = fid->m_Data->GetPixelContainer()->GetBufferPointer();
  OutputPrecisionType* blockPointer = reinterpret_cast< OutputPrecisionType* >( dataPointer );

  int blockCount = 0;
  for( int i = 0; i < header.nblocks; ++i )
    {
    datablockhead blockhead;
    in.read( reinterpret_cast< char* >( &blockhead ), sizeof( datablockhead ) );
    if ( in.fail() || in.eof() )
      {
      break;
      }

    DATABLOCKHEADER_CONVERT_NTOH( &blockhead );

    // allocate block

    bool zero = true;
    for( int j = 0; j < header.ntraces; ++j )
      {
      in.read( reinterpret_cast< char* >( buffer ), sizeof( InputPrecisionType ) * header.np );
      if ( in.fail() || in.bad() || in.eof() )
        {
        break;
        }

      for( int k = 0; k < header.np; ++k, ++blockPointer )
        {
        zero = zero && buffer[ k ] == 0;
        ( *blockPointer ) = static_cast< OutputPrecisionType >( buffer[ k ] );
        ByteSwap( reinterpret_cast< char* >( blockPointer ), sizeof( OutputPrecisionType ) );
        }
      }

    if ( zero )
      {
      break;
      }

    ++blockCount;
    }

  header.nblocks = blockCount;
  size[ 2 ] = blockCount;
  region.SetSize( size );
  fid->m_Data->SetRegions( region );
}

void FIDReader::Read()
{
  m_FID = FID::New();

  datafilehead& header = m_FID->m_Impl->Header;

  std::string path = this->GetFullPath( "fid" );
  std::ifstream in( path.c_str() );
  if ( !in.good() )
    {
    fidThrowMacro( "Could not read from " << path );
    return;
    }

  in.read( reinterpret_cast< char* >( &header ), sizeof( datafilehead ) );

  DATAFILEHEADER_CONVERT_NTOH( &header );

  if ( header.ebytes == 2 )
    {
    ReadFID< short >( in, m_FID );
    }
  else
    {
    ReadFID< float >( in, m_FID );
    }

  m_FID->m_Impl->Procpar.Parse( this->GetFullPath( "procpar" ) );
}

FID::Pointer FIDReader::GetFID()
{
  return m_FID;
}

} // end namespace fid
