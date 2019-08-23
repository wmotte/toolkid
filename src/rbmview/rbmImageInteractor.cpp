#include "rbmView.h"
#include "rbmThreeView.h"

#include "rbmImageInteractor.h"
namespace rbm
{

ImageInteractor::ImageInteractor( ThreeView* view, View* viewer )
: m_View( view ), m_Viewer( viewer )
{

}

void ImageInteractor::SelectCoordinate( double x, double y, double z )
{
	Image::CoordinateType point;
	point[ 0 ] = x;
	point[ 1 ] = y;
	point[ 2 ] = z;
	point[ 3 ] = m_View->GetPlane( 0 )->GetCoordinate()[ 3 ];

	for( int i = 0; i < 3; ++i )
		{
		m_View->GetPlane( i )->SetCoordinate( point );
		}

	m_View->SelectVoxel();
}

void ImageInteractor::SelectVoxel( int x, int y, int z )
{
//	ImagePlane::IndexType voxel;
//	voxel[ 0 ] = x;
//	voxel[ 1 ] = y;
//	voxel[ 2 ] = z;
//
//	for( int i = 0; i < 3; ++i )
//		{
//		m_View->GetPlane( i )->SetCoordinate( m_View->GetPlane( i )->VoxelToCoordinate( voxel ) );
//		}
//
//	m_View->SelectVoxel( x, y, z );
}

void ImageInteractor::SelectVolume( int v )
{
	m_Viewer->SetVolume( v );
}

void ImageInteractor::Plot( Image::VoxelType start, Image::VoxelType end )
{
	m_Viewer->Plot( start, end );
}

void ImageInteractor::TogglePlot()
{
	m_Viewer->TogglePlot();
}

Image::Pointer ImageInteractor::GetImage()
{
  return m_Viewer->m_Image;
}

int ImageInteractor::GetSelectedVolume()
{
  return this->GetImage()->GetVolume();
}

} // end namespace rbm
