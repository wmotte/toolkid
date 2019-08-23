#include "rbmImagePlane.h"
#include "vtkRenderWindow.h"
#include "vtkImageData.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkImageActor.h"
#include "vtkCommand.h"
#include "vnl/vnl_math.h"
#include "rbmImageInteractor.h"
#include "vtkCubeSource.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkSource.h"

namespace rbm
{

ImagePlane::ImagePlane( QWidget *parent ) :
  QWidget( parent ),
  m_Direction( 2 )
{
  ui.setupUi( this );
  ui.m_VTK->setAttribute( Qt::WA_NoSystemBackground );
}

ImagePlane::~ImagePlane()
{

}

void ImagePlane::AddImage( Image::Pointer image, ImageInteractor* imageInteractor )
{
  vtkSmartPointer< vtkImageMapToWindowLevelColors > lut = vtkSmartPointer< vtkImageMapToWindowLevelColors >::New();
  vtkSmartPointer< vtkImageActor > actor = vtkSmartPointer< vtkImageActor >::New();

  if ( m_Images.size() > 0 )
  	{
  	actor->SetOpacity( 0.5 );
  	}

  int extent[ 6 ];
  image->GetImageData()->GetExtent( extent );

  Image::VoxelType voxel;

  if ( m_Images.size() == 1 )
  	{
		for( int i = 0; i < 3; ++i )
			{
			voxel[ i ] = ( extent[ 2 * i ] + extent[ 2 * i + 1 ] ) / 2;
			}
		voxel[ 3 ] = 0;
		m_Coordinate = image->VoxelToCoordinate( voxel );
  	}
  else
  	{
  	voxel = image->CoordinateToVoxel( m_Coordinate );
  	}

  m_Images.push_back( image );
  m_Actors.push_back( actor );
  m_Luts.push_back( lut );
  m_Voxels.push_back( voxel );

  double range[ 2 ];
  image->GetImageData()->GetScalarRange( range );

  double ravg = 0.5 * ( range[ 1 ] + range[ 0 ] );
  double rpro = 0.49 * ( range[ 1 ] - range[ 0 ] );

//	vtkSmartPointer< vtkLookupTable > lut = vtkSmartPointer< vtkLookupTable >::New();
//  lut->SetHueRange( 0, 2.0 / 3.0 );
//  lut->SetValueRange( 1, 1 );
//  lut->SetSaturationRange( 1, 1 );
//	lut->SetTableRange( range[ 0 ], range[ 1 ] );
//	lut->Build();

	lut->SetInput( image->GetImageData() );
//	m_Lut->SetLookupTable( lut );
  lut->SetWindow( rpro * 2 );
  lut->SetLevel( ravg );

  actor->SetInput( lut->GetOutput() );
  actor->InterpolateOff();

	vtkRenderer* renderer = this->GetRenderer();
	vtkCamera* camera = renderer->GetActiveCamera();

	if ( m_Images.size() == 1 )
  	{
		camera->ParallelProjectionOn();
  	}

  renderer->AddActor( actor );

  if ( m_Images.size() == 1 )
  	{
  	this->SetDirection( 2 );
  	}

  if ( !m_InteractorStyle )
  	{
		m_InteractorStyle = vtkSmartPointer< InteractorStyle >::New();
		m_InteractorStyle->SetInteractor( this->GetInteractor() );
		m_InteractorStyle->SetRenderer( renderer );
		m_InteractorStyle->SetImagePlane( this );
	  this->GetInteractor()->SetInteractorStyle( m_InteractorStyle );
  	}

	m_InteractorStyle->SetImageInteractor( imageInteractor );
  m_InteractorStyle->SetImage( image, m_Images.size() - 1 );
}

int ImagePlane::GetDirection() const
{
	return m_Direction;
}

int ImagePlane::GetSlice( int image ) const
{
	return image >= m_Voxels.size() ? 0 : m_Voxels[ image ][ m_Direction ];
}

const Image::VoxelType& ImagePlane::GetVoxel( int image ) const
{
	return m_Voxels[ image ];
//	return image >= m_Voxels.size() ? Image::VoxelType() : m_Voxels[ image ];
}

const Image::CoordinateType& ImagePlane::GetCoordinate() const
{
	return m_Coordinate;
}

void ImagePlane::SetSlice( int slice, int image  )
{
	Image::VoxelType voxel = m_Voxels[ image ];
	voxel[ m_Direction ] = slice;

	this->SetCoordinate( m_Images[ image ]->VoxelToCoordinate( voxel ) );
}

void ImagePlane::SetCoordinate( const Image::CoordinateType& coordinate )
{
	m_Coordinate = coordinate;

	for( int i = 0; i < m_Images.size(); ++i )
		{
		vtkSmartPointer< vtkImageActor > actor = m_Actors[ i ];
		Image::Pointer source = m_Images[ i ];
		vtkSmartPointer< vtkImageData > image = source->GetImageData();
		Image::VoxelType& voxel = m_Voxels[ i ];
		voxel = source->CoordinateToVoxel( coordinate );
		int extent[ 6 ];
		image->GetWholeExtent( extent );

		if ( voxel[ m_Direction ] > extent[ m_Direction * 2 + 1 ] )
			{
			voxel[ m_Direction ] = extent[ m_Direction * 2 + 1 ];
			}
		else if ( voxel[ m_Direction ] < extent[ m_Direction * 2 ] )
			{
			voxel[ m_Direction ] = extent[ m_Direction * 2 ];
			}

		extent[ m_Direction * 2 ] = voxel[ m_Direction ];
		extent[ m_Direction * 2 + 1 ] = voxel[ m_Direction ];

	  actor->SetDisplayExtent( extent[ 0 ], extent[ 1 ], extent[ 2 ], extent[ 3 ], extent[ 4 ], extent[ 5 ] );
		}

	vtkRenderer* renderer = this->GetRenderer();
	vtkCamera* camera = renderer->GetActiveCamera();

	renderer->ResetCameraClippingRange();
	this->GetInteractor()->Render();

//	if ( m_Direction == 2 )
//		{
//		for( int i = 0; i < m_Actors.size(); ++i )
//			{
//			std::cout << "Actor " << i << ":" << std::endl;
//			m_Actors[ i ]->Print( std::cout );
//			}
//		}
}

void ImagePlane::SetDirection( int direction )
{
	m_Direction = direction;

	this->SetCoordinate( m_Coordinate );

	vtkRenderer* renderer = this->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
	vtkCamera* camera = renderer->GetActiveCamera();

  switch( m_Direction )
		{
		case 2:
			{
			camera->SetViewUp( 0, 1, 0 );
			camera->SetPosition( 0, 0, -1 ); // or 0 and focal 1
			camera->SetFocalPoint( 0, 0, 0 );
			break;
			}

		case 1:
			{
			camera->SetViewUp( 1, 0, 0 );
			camera->SetPosition( 0, 1, 0 );
			camera->SetFocalPoint( 0, 0, 0 );
			break;
			}

		case 0:
			{
			camera->SetViewUp( 0, -1, 0 );
			camera->SetPosition( -1, 0, 0 );
			camera->SetFocalPoint( 0, 0, 0 );
			break;
			}
		}

	renderer->ResetCamera();
}

vtkRenderWindow* ImagePlane::GetRenderWindow()
{
  return this->GetWidget()->GetRenderWindow();
}

QVTKInteractor* ImagePlane::GetInteractor()
{
  return this->GetWidget()->GetInteractor();
}

vtkRenderer* ImagePlane::GetRenderer()
{
	return this->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
}

QVTKWidget* ImagePlane::GetWidget()
{
  return ui.m_VTK;
}

Image::CoordinateType ImagePlane::GetEventCoordinate( int image )
{
	vtkRenderWindowInteractor* interactor = this->GetInteractor();
	vtkRenderer* renderer = this->GetRenderer();

	int* pos = interactor->GetEventPosition();
	double world[ 4 ];
	renderer->SetDisplayPoint( pos[0], pos[1], renderer->GetZ( pos[0], pos[1] ) );
	renderer->DisplayToWorld();
	renderer->GetWorldPoint( world );

	Image::CoordinateType point;
	point[ 0 ] = world[ 0 ];
	point[ 1 ] = world[ 1 ];
	point[ 2 ] = world[ 2 ];
	point[ 3 ] = m_Images[ image ]->GetZ();

	double spacing[ 3 ];
	double origin[ 3 ];

	m_Images[ image ]->GetImageData()->GetSpacing( spacing );
  m_Images[ image ]->GetImageData()->GetOrigin( origin );

	point[ m_Direction ] = m_Voxels[ image ][ m_Direction ] * spacing[ m_Direction ] + origin[ m_Direction ];

	return point;
}

float ImagePlane::GetValue( int image )
{
  return m_Images[ image ]->GetImageData()->GetScalarComponentAsFloat( m_Voxels[ image ][ 0 ], m_Voxels[ image ][ 1 ], m_Voxels[ image ][ 2 ], 0 );
}

double ImagePlane::GetWindow( int image ) const
{
  return this->m_Luts[ image ]->GetWindow();
}

double ImagePlane::GetLevel( int image ) const
{
  return this->m_Luts[ image ]->GetLevel();
}

void ImagePlane::SetWindowLevel( double window, double level, int image  )
{
  this->m_Luts[ image ]->SetLevel( level );
  this->m_Luts[ image ]->SetWindow( window );
  this->GetInteractor()->Render();
}

} // end namespace rbm
