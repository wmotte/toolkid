PROJECT( rbmview )

if( VTK_FOUND )
if( QT_FOUND )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )

QT4_WRAP_UI( rbmView_UIS_H
	rbmView.ui
	rbmThreeView.ui
	rbmImagePlane.ui
)

QT4_WRAP_CPP( rbmView_MOC_SRCS
	rbmView.h
	rbmThreeView.h
	rbmImagePlane.h
	rbmTest.h
)

ADD_EXECUTABLE( rbmview
	WIN32
	MACOSX_BUNDLE
	${rbmView_UIS_H}
	${rbmView_MOC_SRCS}
	main.cpp
	rbmThreeView.cpp
	rbmImagePlane.cpp
	rbmView.cpp
	rbmImageInteractor.h
	rbmImageInteractor.cpp
	rbmPlotCurve.h
	rbmPlotCurve.cpp
	rbmPlotData.h
	rbmPlotData.cpp
	rbmPlotter.h
	rbmPlotter.cpp
	rbmImage.h
	rbmImage.cpp
	rbmInteractorStyle.h
	rbmInteractorStyle.cpp
)

TARGET_LINK_LIBRARIES( rbmview
	QVTK
	qwt-qt4
	${QT_LIBRARIES}
	ITKCommon
	ITKIO
	ITKBasicFilters
	vtkRendering
	vtkGraphics
	vtkIO
	vtkCommon
)

endif( QT_FOUND )
endif( VTK_FOUND )
