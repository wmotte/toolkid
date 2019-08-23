#include "fidFID.h"
#include "fidFIDImpl.h"

namespace fid
{

FID::FID() : m_Impl( new FIDImpl, true )
{

}

FID::~FID()
{

}

FID::DataPointer FID::GetData()
{
	return m_Data;
}

FID::DataConstPointer FID::GetData() const
{
	return m_Data.GetPointer();
}

FID::BlockType FID::GetBlock( int block )
{
	DataType::IndexType index;
	index[ 0 ] = 0;
	index[ 1 ] = 0;
	index[ 2 ] = block;

	ComplexType& pixel = m_Data->GetPixel( index );

  return BlockType( this->GetNumberOfTraces(), this->GetNumberOfPoints(), &pixel );
}

int FID::GetNumberOfBlocks() const
{
  return m_Impl->Header.nblocks;
}

int FID::GetNumberOfTraces() const
{
  return m_Impl->Header.ntraces;
}

int FID::GetNumberOfPoints() const
{
  return m_Impl->Header.np / 2;
}

const Procparser& FID::GetProcpar() const
{
  return m_Impl->Procpar;
}

} // end namespace fid
