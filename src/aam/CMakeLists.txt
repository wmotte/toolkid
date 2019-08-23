cmake_minimum_required( VERSION 2.6 )
 
project( aam )
    
find_package( VTK REQUIRED )  
include( ${VTK_USE_FILE} )
   
# procrustes test ...
#add_executable( ProcrustesAlignmentFilter ProcrustesAlignmentFilter.cxx )
#target_link_libraries( ProcrustesAlignmentFilter vtkHybrid )

# active shape model test ...
add_executable( aam_test itkImagePCAShapeModelEstimatorTest.cxx )
target_link_libraries( aam_test ITKCommon ITKIO ITKStatistics ITKBasicFilters )

# read/write vtk ...
#add_executable( aam_vtk itkVTKImageIOTest.cxx )
#target_link_libraries( aam_vtk ITKCommon ITKIO ITKStatistics ITKBasicFilters )
