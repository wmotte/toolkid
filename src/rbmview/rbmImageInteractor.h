#ifndef __rbmImageInteractor_h__
#define __rbmImageInteractor_h__
#include "itkIndex.h"
#include "rbmImage.h"

namespace rbm
{

class ThreeView;
class View;

class ImageInteractor
{
public:
	ImageInteractor( ThreeView* view, View* viewer );

	void SelectCoordinate( double x, double y, double z );
	void SelectVoxel( int x, int y, int z );
	void SelectVolume( int v );
	void Plot( Image::VoxelType start, Image::VoxelType end );
	void TogglePlot();
	Image::Pointer GetImage();
	int GetSelectedVolume();

protected:
	ThreeView* m_View;
	View* m_Viewer;
	bool m_Plot;
};

} // end namespace rbm

#endif /*__rbmImageInteractor_h__*/
