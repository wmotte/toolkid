/****************************************************************************
**
** Trolltech hereby grants a license to use the Qt/Eclipse Integration
** plug-in (the software contained herein), in binary form, solely for the
** purpose of creating code to be used with Trolltech's Qt software.
**
** Qt Designer is licensed under the terms of the GNU General Public
** License versions 2.0 and 3.0 ("GPL License"). Trolltech offers users the
** right to use certain no GPL licensed software under the terms of its GPL
** Exception version 1.2 (http://trolltech.com/products/qt/gplexception).
**
** THIS SOFTWARE IS PROVIDED BY TROLLTECH AND ITS CONTRIBUTORS (IF ANY) "AS
** IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
** TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
** PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
** OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
** EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
** PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** Since we now have the GPL exception I think that the "special exception
** is no longer needed. The license text proposed above (other than the
** special exception portion of it) is the BSD license and we have added
** the BSD license as a permissible license under the exception.
**
****************************************************************************/

#ifndef RBMIMAGEPLANE_H
#define RBMIMAGEPLANE_H

#include <QtGui/QWidget>
#include "QVTKWidget.h"
#include "ui_rbmImagePlane.h"
#include "vtkSmartPointer.h"
#include "rbmImage.h"
#include "rbmInteractorStyle.h"
#include <vector>

class vtkRenderWindow;
class QVTKInteractor;
class vtkImageData;
class vtkImageActor;
class vtkImageMapToWindowLevelColors;

namespace rbm
{

class ImageInteractor;

class ImagePlane : public QWidget
{
    Q_OBJECT

public:
    ImagePlane( QWidget *parent = 0 );
    ~ImagePlane();

    vtkRenderWindow* GetRenderWindow();
    vtkRenderer* GetRenderer();
    QVTKInteractor* GetInteractor();
    QVTKWidget* GetWidget();

    void AddImage( Image::Pointer image, ImageInteractor* imageInteractor );
    void SetDirection( int direction );
    void SetCoordinate( const Image::CoordinateType& coordinate );
    void SetSlice( int slice, int image = 0 );
    const Image::VoxelType& GetVoxel( int image = 0 ) const;
    const Image::CoordinateType& GetCoordinate() const;
    int GetSlice( int image = 0) const;
    int GetDirection() const;
    double GetWindow( int image = 0 ) const;
    double GetLevel( int image = 0 ) const;
    void SetWindowLevel( double window, double level, int image = 0 );

    Image::CoordinateType GetEventCoordinate( int image = 0 );
    float GetValue( int image = 0 );

protected:
    std::vector< Image::Pointer > m_Images;
    std::vector< Image::VoxelType > m_Voxels;
    std::vector< vtkSmartPointer< vtkImageActor > > m_Actors;
    std::vector< vtkSmartPointer< vtkImageMapToWindowLevelColors > > m_Luts;
    vtkSmartPointer< InteractorStyle > m_InteractorStyle;
    int m_Direction;
    Image::CoordinateType m_Coordinate;

private:
    Ui::ImagePlaneClass ui;
};

} // end namespace rbm

#endif // RBMIMAGEPLANE_H
