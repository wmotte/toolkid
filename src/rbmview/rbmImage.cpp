#include "vtkImageData.h"
#include "rbmImage.h"
#include <vector>

namespace rbm
{

Image::Image()
  : m_Matrix( 0, 0, 0 )
{
}

void Image::SetSeries( SeriesType::Pointer series, int volume )
{
	this->Initialize( series, volume );
}

void Image::Initialize( SeriesType::Pointer series, int volume )
{
	m_Series = series;

	if ( !m_Series )
		{
		m_Series = SeriesType::New();
		m_Size.Fill( 1 );
		m_Region.SetSize( m_Size );
		m_Series->SetRegions( m_Region );
		m_Series->Allocate();
		m_Series->FillBuffer( 0 );
		}

//	m_Series->Update();

	m_Region = m_Series->GetLargestPossibleRegion();
	m_Size = m_Region.GetSize();
	m_Spacing = m_Series->GetSpacing();
	m_Origin = m_Series->GetOrigin();

	m_ExtractionRegion = m_Region;
	m_ExtractionSize = m_Size;
	m_ExtractionIndex = m_Region.GetIndex();

	m_Volume = volume;
	m_Volumes = m_Size[ 3 ];
	m_Voxels = m_Size[ 0 ] * m_Size[ 1 ] * m_Size[ 2 ];

	m_ExtractionIndex[ 3 ] = m_Volume;
	m_ExtractionSize[ 3 ] = 0;
	m_ExtractionRegion.SetSize( m_ExtractionSize );
	m_ExtractionRegion.SetIndex( m_ExtractionIndex );

	m_Z = this->VoxelToCoordinate( m_ExtractionIndex )[ 3 ];

	m_Extract = ExtractType::New();
	m_Image = m_Extract->GetOutput();

	m_Extract->SetInput( m_Series );

	m_Filter = FilterType::New();
	m_Filter->SetInput( m_Extract->GetOutput() );
	m_ImageData = m_Filter->GetOutput();

	m_Extract->SetExtractionRegion( m_ExtractionRegion );

	Update();
}

void Image::Update()
{
	m_ExtractionIndex[ 3 ] = m_Volume;
	m_ExtractionRegion.SetIndex( m_ExtractionIndex );
	m_Extract->SetExtractionRegion( m_ExtractionRegion );
	m_Z = this->VoxelToCoordinate( m_ExtractionIndex )[ 3 ];
  m_Extract->Update();
	m_Filter->Update();
////	m_Extract->UpdateLargestPossibleRegion();
//	m_Filter->UpdateLargestPossibleRegion();
}

Image::ImageType::Pointer Image::GetImage() const
{
	return m_Image;
}

vtkSmartPointer< vtkImageData > Image::GetImageData() const
{
	return m_ImageData;
}

Image::MatrixType Image::GetImageMatrix() const
{
  return MatrixType( m_Volumes, m_Voxels, m_Series->GetPixelContainer()->GetBufferPointer() );
}

Image::SeriesType::Pointer Image::GetSeries() const
{
	return m_Series;
}

Image::SpacingType Image::GetSpacing() const
{
	return m_Spacing;
}

Image::CoordinateType Image::GetOrigin() const
{
	return m_Origin;
}

int Image::GetNumberOfVolumes() const
{
	return m_Volumes;
}

int Image::GetNumberOfVoxels() const
{
	return m_Voxels;
}

Image::SizeType Image::GetSize() const
{
	return m_Size;
}

void Image::SetVolume( int volume )
{
  if ( volume >= this->GetNumberOfVolumes() )
    {
    volume = 0;
    }
  else if ( volume < 0 )
    {
    volume = this->GetNumberOfVolumes() - 1 ;
    }

	m_Volume = volume;
	Update();
}

int Image::GetVolume() const
{
	return m_Volume;
}

double Image::GetZ() const
{
	return m_Z;
}

Image::VoxelType Image::CoordinateToVoxel( const CoordinateType& coordinate ) const
{
	VoxelType voxel;

	for( int i = 0; i < 4; ++i )
		{
		voxel[ i ] = vcl_ceil( ( coordinate[ i ] - m_Origin [ i ] - m_Spacing[ i ] * 0.5 ) / m_Spacing[ i ] );
		}

	return voxel;
}

Image::CoordinateType Image::VoxelToCoordinate( const VoxelType& voxel ) const
{
	CoordinateType coordinate;

	for( int i = 0; i < 4; ++i )
		{
		coordinate[ i ] = static_cast< double >( voxel[ i ] ) * m_Spacing[ i ] + m_Origin[ i ];
		}

	return coordinate;
}

int Image::GetProfile( VoxelType start, VoxelType end, VectorType& mean, VectorType& sd ) const
{
	mean = VectorType( m_Volumes );
	sd = VectorType( m_Volumes );

	mean.fill( 0 );
	sd.fill( 0 );
	int n = 0;

	std::vector< int > indices;

	MatrixType matrix( m_Volumes, m_Voxels, m_Series->GetPixelContainer()->GetBufferPointer() );

	for( int i = 0; i < 3; ++i )
		{
		if ( start[ i ] >= m_Size[ i ] )
			{
			start[ i ] = m_Size[ i ] - 1;
			}

		if ( end[ i ] >= m_Size[ i ] )
			{
			end[ i ] = m_Size[ i ] - 1;
			}

		if ( start[ i ] < 0  )
			{
			start[ i ] = 0;
			}

		if ( end[ i ] < 0  )
			{
			end[ i ] = 0;
			}

		if ( start[ i ] > end[ i ] )
			{
			int t = end[ i ];
			end[ i ] = start[ i ];
			start[ i ] = t;
			}
		}

	for( int x = start[ 0 ]; x <= end[ 0 ]; ++x )
    {
    for( int y = start[ 1 ]; y <= end[ 1 ]; ++y )
      {
      for( int z = start[ 2 ]; z <= end[ 2 ]; ++z )
        {
				int index = z * m_Size[ 1 ] * m_Size[ 0 ] + y * m_Size[ 0 ] + x;
				indices.push_back( index );

				mean += matrix.get_column( index );
				++n;
        }
      }
    }

  mean /= static_cast< PixelType >( n );

  for( std::vector< int >::iterator i = indices.begin(); i != indices.end(); ++i )
  	{
  	VectorType profile = matrix.get_column( *i ) - mean;

  	for( int j = 0; j < m_Volumes; ++j )
  		{
  		sd[ j ] += profile[ j ] * profile[ j ];
  		}
  	}

  for( int i = 0; i < m_Volumes; ++i )
  	{
  	sd[ i ] = vcl_sqrt( sd[ i ] / static_cast< PixelType >( n - 1 ) );
  	}

  return n;
}

} // end namespace rbm
