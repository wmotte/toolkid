#ifndef RBMTHREEVIEW_H
#define RBMTHREEVIEW_H

#include <QtGui/QWidget>

#include "ui_rbmThreeView.h"
#include "vtkSmartPointer.h"
#include "vtkRenderer.h"
#include "rbmImagePlane.h"
#include "rbmImage.h"

#include <QApplication>


namespace rbm
{

class ImageInteractor;

class ThreeView : public QWidget
{
  Q_OBJECT

public:
  ThreeView( QWidget *parent = 0 );
  ~ThreeView();

  vtkRenderer* GetRenderer( int index );
  ImagePlane* GetPlane( int index );

  void SelectVoxel();

  void AddImage( Image::Pointer image, ImageInteractor* interactor, int v );
  void SelectImage( int index );
  int GetSelectedImage() const;
  Image::Pointer GetImage( int index ) const;
  void Render();

public slots:
  void UpdateValueX( double x );
  void UpdateValueY( double y );
  void UpdateValueZ( double z );
  void UpdateValueT( int t );
  void UpdateWindow( double window );
  void UpdateLevel( double level );
  void ViewGraph( void );

protected:
  vtkSmartPointer< vtkRenderer > m_Renderer[ 3 ];
  ImagePlane* m_Widget[ 3 ];
  ImageInteractor* m_ImageInteractor;
  std::vector< Image::Pointer > m_Images;
  int m_ActiveImage;

private:
  Ui::ThreeViewClass ui;
};

} // end namespace rbm

#endif // RBMTHREEVIEW_H
