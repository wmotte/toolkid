#ifndef __rbmInteractorStyle_h__
#define __rbmInteractorStyle_h__
#include "rbmImage.h"
#include "vtkInteractorStyleImage.h"

class vtkRenderWindowInteractor;
class vtkRenderer;
class vtkActor;
class vtkPolyDataMapper;
class vtkCubeSource;

namespace rbm
{

class ImagePlane;
class ImageInteractor;

class InteractorStyle : public vtkInteractorStyleImage
{
public:
	static InteractorStyle* New();

	void SetInteractor( vtkSmartPointer< vtkRenderWindowInteractor > interactor );
	void SetImagePlane( ImagePlane* plane );
	void SetRenderer( vtkSmartPointer< vtkRenderer > renderer );
	void SetImage( Image::Pointer image, int index );
	void SetImageInteractor( ImageInteractor* interactor );

	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnMouseMove();
  virtual void OnChar();
  virtual void OnKeyPress();

protected:
	InteractorStyle();

	void ShowBox();
	void HideBox();
	void SelectVoxel();
	void Plot( Image::CoordinateType start, Image::CoordinateType end );

	vtkSmartPointer< vtkActor > m_BoxActor;
	vtkSmartPointer< vtkPolyDataMapper > m_BoxMapper;
	vtkSmartPointer< vtkCubeSource > m_BoxSource;

	bool m_Select;
	bool m_LeftDown;
	bool m_BoxVisible;
	int m_ActiveImage;
	Image::CoordinateType m_Start;
	Image::CoordinateType m_End;

	Image::Pointer m_Image;
	vtkSmartPointer< vtkRenderWindowInteractor > m_Interactor;
	vtkSmartPointer< vtkRenderer > m_Renderer;
	ImagePlane* m_Plane;
  ImageInteractor* m_ImageInteractor;

  double m_Spacing[ 3 ];
	double m_Origin[ 3 ];
	int m_Extent[ 6 ];
};

} // end namespace rbm

#endif /*__rbmInteractorStyle_h__*/
