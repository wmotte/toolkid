#include "rbmInteractorStyle.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkCubeSource.h"
#include "rbmImagePlane.h"
#include "rbmImageInteractor.h"
#include "vtkObjectFactory.h"
#include "vtkProperty.h"

namespace rbm
{

vtkStandardNewMacro( InteractorStyle );

InteractorStyle::InteractorStyle()
{
	m_Select = false;
	m_LeftDown = false;
	m_BoxActor = vtkSmartPointer< vtkActor >::New();
	m_BoxMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
	m_BoxSource = vtkSmartPointer< vtkCubeSource >::New();
	m_BoxVisible = false;

	m_BoxMapper->SetInputConnection( m_BoxSource->GetOutputPort() );
	m_BoxActor->SetMapper( m_BoxMapper );
//		m_BoxActor->GetProperty()->SetRepresentationToWireframe();
	m_BoxActor->GetProperty()->SetOpacity( 0.5 );
	m_BoxActor->GetProperty()->SetColor( 1, 1, 0 );
}

void InteractorStyle::OnLeftButtonDown()
{
	m_Select = false;
	m_LeftDown = true;
	m_Start = m_Image->VoxelToCoordinate( m_Image->CoordinateToVoxel( m_Plane->GetEventCoordinate( m_ActiveImage ) ) );
	HideBox();

	vtkInteractorStyleImage::OnLeftButtonDown();
}

void InteractorStyle::OnLeftButtonUp()
{
	if ( !m_Select )
		{
		SelectVoxel();
		}
	else
		{
		m_End = m_Plane->GetEventCoordinate( m_ActiveImage );
//		std::cout << "Selected " << m_Start << " - " << end << std::endl;

//		HideBox();

		Plot( m_Start, m_End );
		}

	m_LeftDown = false;
	m_Select = false;

	vtkInteractorStyleImage::OnLeftButtonUp();
}

void InteractorStyle::ShowBox()
{
	if ( !m_BoxVisible )
		{
		m_Renderer->AddActor( m_BoxActor );
		m_BoxVisible = true;
		}
}

void InteractorStyle::HideBox()
{
	if ( m_BoxVisible )
		{
		m_Renderer->RemoveActor( m_BoxActor );
		m_BoxVisible = false;
		m_Renderer->ResetCameraClippingRange();
		m_Interactor->Render();
		}
}

void InteractorStyle::OnMouseMove()
{
	if ( !m_LeftDown && !m_Select )
		{
		vtkInteractorStyleImage::OnMouseMove();
		return;
		}

	m_Select = true;

	if ( !m_BoxVisible )
		{
		ShowBox();
		}

	Image::CoordinateType start = m_Start;
	Image::CoordinateType end = m_Image->VoxelToCoordinate( m_Image->CoordinateToVoxel( m_Plane->GetEventCoordinate( m_ActiveImage ) ) );//m_Plane->GetEventCoordinate();

	for( int i = 0; i < 3; ++i )
		{
		if ( start[ i ] < end[ i ] )
			{
			start[ i ] -= 0.5 * m_Spacing[ i ];
			end[ i ] += 0.5 * m_Spacing[ i ];
			}
		else
			{
			start[ i ] += 0.5 * m_Spacing[ i ];
			end[ i ] -= 0.5 * m_Spacing[ i ];
			}
		}

	double bounds[ 6 ];

	for( int i = 0; i < 3; ++i )
		{
		if ( start[ i ] == end[ i ] )
			{
			start[ i ] -= 0.1;
			end[ i ] += 0.1;
			}

		if ( start[ i ] > end[ i ] )
			{
			bounds[ i * 2 ] = end[ i ];
			bounds[ i * 2 + 1 ] = start[ i ];
			}
		else
			{
			bounds[ i * 2 ] = start[ i ];
			bounds[ i * 2 + 1 ] = end[ i ];
			}
		}

	m_BoxSource->SetBounds( bounds );
	m_Renderer->ResetCameraClippingRange();
	m_Interactor->Render();

	vtkInteractorStyleImage::OnMouseMove();
}

void InteractorStyle::OnChar()
{
  vtkRenderWindowInteractor* rwi = this->Interactor;
  switch( rwi->GetKeyCode() )
    {
    case '+':
    case '=':
      {
      this->m_ImageInteractor->SelectVolume( this->m_ImageInteractor->GetSelectedVolume() + 1 );
      break;
      }

    case '-':
    case '_':
      {
      this->m_ImageInteractor->SelectVolume( this->m_ImageInteractor->GetSelectedVolume() - 1 );
      break;
      }
    }
}

void InteractorStyle::OnKeyPress()
{
  vtkRenderWindowInteractor* rwi = this->Interactor;
  std::string sym = rwi->GetKeySym();

  if ( sym == "Up" )
    {
    m_Start[ 1 ] += m_Spacing[ 1 ];
    m_End[ 1 ] += m_Spacing[ 1 ];
    }
  else if ( sym == "Down" )
    {
    m_Start[ 1 ] -= m_Spacing[ 1 ];
    m_End[ 1 ] -= m_Spacing[ 1 ];
    }
  else if ( sym == "Left" )
    {
    m_Start[ 0 ] += m_Spacing[ 0 ];
    m_End[ 0 ] += m_Spacing[ 0 ];
    }
  else if ( sym == "Right" )
    {
    m_Start[ 0 ] -= m_Spacing[ 0 ];
    m_End[ 0 ] -= m_Spacing[ 0 ];
    }
  else
    {
    return;
    }

  Plot( m_Start, m_End );

  Image::CoordinateType start = m_Start;
  Image::CoordinateType end = m_End;

  for( int i = 0; i < 3; ++i )
    {
    if ( start[ i ] < end[ i ] )
      {
      start[ i ] -= 0.5 * m_Spacing[ i ];
      end[ i ] += 0.5 * m_Spacing[ i ];
      }
    else
      {
      start[ i ] += 0.5 * m_Spacing[ i ];
      end[ i ] -= 0.5 * m_Spacing[ i ];
      }
    }

  double bounds[ 6 ];

  for( int i = 0; i < 3; ++i )
    {
    if ( start[ i ] == end[ i ] )
      {
      start[ i ] -= 0.1;
      end[ i ] += 0.1;
      }

    if ( start[ i ] > end[ i ] )
      {
      bounds[ i * 2 ] = end[ i ];
      bounds[ i * 2 + 1 ] = start[ i ];
      }
    else
      {
      bounds[ i * 2 ] = start[ i ];
      bounds[ i * 2 + 1 ] = end[ i ];
      }
    }

  m_BoxSource->SetBounds( bounds );
  m_Renderer->ResetCameraClippingRange();
  m_Interactor->Render();
}

void InteractorStyle::SelectVoxel()
{
	Image::CoordinateType point = m_Plane->GetEventCoordinate( m_ActiveImage );
//	std::cout << "Point: " << point << std::endl;
	m_ImageInteractor->SelectCoordinate( point[ 0 ], point[ 1 ], point[ 2 ] );
	Plot( point, point );
}

void InteractorStyle::Plot( Image::CoordinateType start, Image::CoordinateType end )
{
	m_ImageInteractor->Plot( m_Image->CoordinateToVoxel( start ), m_Image->CoordinateToVoxel( end ) );
}

void InteractorStyle::SetInteractor( vtkSmartPointer< vtkRenderWindowInteractor > interactor )
{
	m_Interactor = interactor;
}

void InteractorStyle::SetImagePlane( ImagePlane* plane )
{
	m_Plane = plane;
}

void InteractorStyle::SetRenderer( vtkSmartPointer< vtkRenderer > renderer )
{
	m_Renderer = renderer;
}

void InteractorStyle::SetImage( Image::Pointer image, int index )
{
	m_Image = image;
	m_ActiveImage = index;
	m_Image->GetImageData()->GetOrigin( m_Origin );
	m_Image->GetImageData()->GetSpacing( m_Spacing );
	m_Image->GetImageData()->GetExtent( m_Extent );
}

void InteractorStyle::SetImageInteractor( ImageInteractor* interactor )
{
	m_ImageInteractor = interactor;
}

} // end namespace rbm
