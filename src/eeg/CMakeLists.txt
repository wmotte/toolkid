cmake_minimum_required( VERSION 2.6 )
project( eeg )

include_directories( 
	../cmdparser
	${Boost_INCLUDE_DIRS}
)

add_executable( eeg_txt2nii txt2nii.cpp )
target_link_libraries( eeg_txt2nii cmdparser ${Boost_LIBRARIES} ITKCommon ITKIO ITKAlgorithms ITKBasicFilters)

add_executable( eeg_extractSeries extractSeries.cpp )
target_link_libraries( eeg_extractSeries cmdparser ITKCommon ITKIO ITKAlgorithms ITKBasicFilters)
