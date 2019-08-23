#include "rbmView.h"
#include "itkImageFileReader.h"
#include "rbmImageInteractor.h"
#include "rbmPlotter.h"
#include "qwt_plot_picker.h"

namespace rbm
{

View::View( QWidget *parent )
    : QMainWindow( parent )
{
	ui.setupUi( this );
	m_Plotter = 0;
}

View::~View()
{

}

void View::SetVolume( int v )
{
	m_Image->SetVolume( v );

	for( int i = 0; i < 3; ++i )
		{
		Image::CoordinateType coordinate = ui.widget->GetPlane( i )->GetCoordinate();
		coordinate[ 3 ] = m_Image->GetZ();
		ui.widget->GetPlane( i )->SetCoordinate( coordinate );
		}

	ui.widget->Render();
}

void View::LoadOverlay( const QString& filename )
{
//  typedef itk::ImageFileReader< ImageType > ReaderType;
//  ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName( filename.toStdString().c_str() );
//  reader->Update();
//
//  ImageType::Pointer image = reader->GetOutput();
//  reader = 0;
//
//  FilterType::Pointer filter = FilterType::New();
//  filter->SetInput( image );
//  filter->Update();
//
//  ui.widget->GetPlane( 0 )->AddOverlay( filter->GetOutput() );
//  ui.widget->GetPlane( 1 )->AddOverlay( filter->GetOutput() );
//  ui.widget->GetPlane( 2 )->AddOverlay( filter->GetOutput() );
}

void View::LoadImage( const QString& filename )
{
  typedef itk::ImageFileReader< Image::SeriesType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename.toStdString().c_str() );
  reader->Update();

  m_Image = Image::New();
  m_Image->SetSeries( reader->GetOutput() );

  ImageInteractor* imageInteractor = new ImageInteractor( ui.widget, this );
  ui.widget->AddImage( m_Image, imageInteractor, m_Image->GetNumberOfVolumes() );
}

void View::TogglePlot()
{
	if ( m_Plotter )
	  {
	  if ( m_Plotter->isVisible() )
	    {
	    m_Plotter->hide();
	    }
	  else
	    {
	    m_Plotter->show();
	    }
	  }
}

void View::Plot( Image::VoxelType start, Image::VoxelType end )
{
	Image::VectorType profile;
	Image::VectorType sd;

	int n = m_Image->GetProfile( start, end, profile, sd );

  if ( !m_Plotter )
		{
		m_Plotter = new Plotter( profile, sd );
		m_Plotter->resize( 600, 400 );

		QwtPlotPicker* picker = new QwtPlotPicker( m_Plotter->canvas() );
		QObject::connect( picker, SIGNAL( selected( QwtDoublePoint& ) ), this, SLOT( PlotPointSelected( QwtDoublePoint& ) ) );
		}
	else
		{
		m_Plotter->set( profile, sd );
		}
}

void View::PlotPointSelected( QwtDoublePoint& point )
{
  std::cout << point.x() << ", " << point.y() << std::endl;
}

} // end namespace rbm
