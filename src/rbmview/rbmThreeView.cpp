#include "rbmThreeView.h"
#include "vtkConeSource.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "rbmImagePlane.h"
#include "rbmImageInteractor.h"
#include "vtkImageData.h"

namespace rbm
{

ThreeView::ThreeView( QWidget *parent )
    : QWidget( parent )
{
	ui.setupUi( this );

	m_Widget[ 0 ] = this->ui.m_Widget1;
  m_Widget[ 1 ] = this->ui.m_Widget2;
  m_Widget[ 2 ] = this->ui.m_Widget3;

  for( int i = 0; i < 3; ++i )
    {
    vtkSmartPointer< vtkRenderer > renderer = vtkSmartPointer< vtkRenderer >::New();
    m_Widget[ i ]->GetRenderWindow()->AddRenderer( renderer );
    m_Renderer[ i ] = renderer;
    }
}

ThreeView::~ThreeView()
{

}

void ThreeView::Render()
{
	m_Widget[ 0 ]->GetInteractor()->Render();
	m_Widget[ 1 ]->GetInteractor()->Render();
	m_Widget[ 2 ]->GetInteractor()->Render();
}

void ThreeView::SelectVoxel()
{
	Image::CoordinateType coordinate = this->GetPlane( 0 )->GetCoordinate();
	Image::VoxelType voxel = this->GetPlane( 0 )->GetVoxel( m_ActiveImage );

	ui.m_SpinX->setValue( coordinate[ 0 ] );
	ui.m_SpinY->setValue( coordinate[ 1 ] );
	ui.m_SpinZ->setValue( coordinate[ 2 ] );

	ui.m_LabelX->setText( QString::number( static_cast< double >( voxel [ 0 ] ), 'f', 2 ) );
  ui.m_LabelY->setText( QString::number( static_cast< double >( voxel [ 1 ] ), 'f', 2 ) );
  ui.m_LabelZ->setText( QString::number( static_cast< double >( voxel [ 2 ] ), 'f', 2 ) );

  float intensity = m_Images[ m_ActiveImage ]->GetSeries()->GetPixel( voxel );

  ui.m_LabelIntensity->setText( QString::number( intensity, 'f', 6 ) );

  ui.m_SpinT->setValue( m_Images[ m_ActiveImage ]->GetVolume() );
}

void ThreeView::SelectImage( int index )
{
	m_ActiveImage = index;

	Image::Pointer image = m_Images[ index ];

	int extent[ 6 ];
	image->GetImageData()->GetWholeExtent( extent );

	this->ui.m_SpinT->setMaximum( image->GetNumberOfVolumes() );
	this->ui.m_SpinX->setMaximum( extent[ 1 ] );
	this->ui.m_SpinY->setMaximum( extent[ 3 ] );
	this->ui.m_SpinZ->setMaximum( extent[ 5 ] );

	this->ui.m_SpinWindow->setValue( this->GetPlane( 0 )->GetWindow( m_ActiveImage ) );
	this->ui.m_SpinLevel->setValue( this->GetPlane( 0 )->GetLevel( m_ActiveImage ) );

	double range[ 2 ];
	image->GetImageData()->GetScalarRange( range );

	this->ui.m_SpinLevel->setMinimum( range[ 0 ] - range[ 1 ] );
	this->ui.m_SpinLevel->setMaximum( range[ 1 ] + range[ 1 ] );

	this->ui.m_SpinWindow->setMinimum( range[ 0 ] - range[ 1 ] );
	this->ui.m_SpinWindow->setMaximum( range[ 1 ] + range[ 1 ] );

	this->ui.m_SpinLevel->setSingleStep( ( range[ 1 ] - range[ 0 ] ) / 128. );
	this->ui.m_SpinWindow->setSingleStep( ( range[ 1 ] - range[ 0 ] ) / 128. );
}

int ThreeView::GetSelectedImage() const
{
	return m_ActiveImage;
}

Image::Pointer ThreeView::GetImage( int index ) const
{
	return m_Images[ index ];
}

vtkRenderer* ThreeView::GetRenderer( int index )
{
  return m_Renderer[ index ];
}

ImagePlane* ThreeView::GetPlane( int index )
{
  return m_Widget[ index ];
}

void ThreeView::UpdateValueX( double x )
{
  Image::CoordinateType point = this->GetPlane( 0 )->GetCoordinate();
  point[ 0 ] = x;
  this->m_ImageInteractor->SelectCoordinate( point[ 0 ], point[ 1 ], point[ 2 ] );
}

void ThreeView::UpdateValueY( double y )
{
  Image::CoordinateType point = this->GetPlane( 0 )->GetCoordinate();
  point[ 1 ] = y;
  this->m_ImageInteractor->SelectCoordinate( point[ 0 ], point[ 1 ], point[ 2 ] );
}

void ThreeView::UpdateValueZ( double z )
{
  Image::CoordinateType point = this->GetPlane( 0 )->GetCoordinate();
  point[ 2 ] = z;
  this->m_ImageInteractor->SelectCoordinate( point[ 0 ], point[ 1 ], point[ 2 ] );
}

void ThreeView::UpdateValueT( int t )
{
  this->m_ImageInteractor->SelectVolume( t );
}

void ThreeView::UpdateWindow( double window  )
{
  for( int i = 0; i < 3; ++i )
    {
    this->GetPlane( i )->SetWindowLevel( window, this->GetPlane( i )->GetLevel( m_ActiveImage ), m_ActiveImage );
    }
}

void ThreeView::UpdateLevel( double level )
{
  for( int i = 0; i < 3; ++i )
    {
    this->GetPlane( i )->SetWindowLevel( this->GetPlane( i )->GetWindow( m_ActiveImage ), level, m_ActiveImage );
    }
}

void ThreeView::ViewGraph( void )
{
	Image::VoxelType start, end;
	start = this->GetPlane( 0 )->GetVoxel();
	end = start;

//	this->m_ImageInteractor->Plot( start, end );
	this->m_ImageInteractor->TogglePlot();

//	std::cout << "Nice" << std::endl;
//	Plot* plot = new Plot;
//  plot->resize(600,400);
//  plot->show();
}

void ThreeView::AddImage( Image::Pointer image, ImageInteractor* interactor, int v )
{
  bool first = m_Images.size() == 0;

  if ( first )
  	{
  	for( int i = 0; i < 3; ++i)
  		{
  		this->GetPlane( i )->SetDirection( i );
  		}
  	}

	m_ImageInteractor = interactor;

	int extent[ 6 ];
	image->GetImageData()->GetExtent( extent );

	for( int i = 0; i < 3; ++i )
		{
		this->GetPlane( i )->AddImage( image, interactor );

		if ( first )
			{
			this->GetPlane( i )->SetDirection( i );
			this->GetPlane( i )->SetCoordinate( this->GetPlane( i )->GetCoordinate() );
			}
		}

	if ( first )
		{
//		QObject::connect( ui.m_SpinX, SIGNAL( valueChanged( double ) ), this, SLOT( UpdateValueX( double ) ) );
//		QObject::connect( ui.m_SpinY, SIGNAL( valueChanged( double ) ), this, SLOT( UpdateValueY( double ) ) );
//		QObject::connect( ui.m_SpinZ, SIGNAL( valueChanged( double ) ), this, SLOT( UpdateValueZ( double ) ) );

		QObject::connect( ui.pushButton, SIGNAL( clicked( void ) ), this, SLOT( ViewGraph( void ) ) );

		QObject::connect( ui.m_SpinT, SIGNAL( valueChanged( int ) ), this, SLOT( UpdateValueT( int ) ) );

		QObject::connect( ui.m_SpinWindow, SIGNAL( valueChanged( double ) ), this, SLOT( UpdateWindow( double ) ) );
		QObject::connect( ui.m_SpinLevel, SIGNAL( valueChanged( double ) ), this, SLOT( UpdateLevel( double ) ) );
		}

	m_Images.push_back( image );
	this->SelectImage( m_Images.size() - 1 );
}

} // end namespace rbm
