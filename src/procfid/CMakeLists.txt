project( procfid )

include_directories(
.
../cmdparser
../procparser 
../mrfit
../dtifit
)

add_library(
	FID
	fidFID.h
	fidFID.cpp
	fidFIDReader.h
	fidFIDReader.cpp
	fidFourier.h
	fidFourier.cpp
	fidEPITools.h
	fidEPITools.cpp
	fidCommon.h
	fidCommon.cpp
)

add_executable(procmsepilcc procmsepilcc.cpp)
target_link_libraries(procmsepilcc FID procparser cmdparser dtifit mrfit ITKCommon ITKIO ITKNumerics)

add_executable(procmsepilc procmsepilc.cpp)
target_link_libraries(procmsepilc FID procparser cmdparser dtifit ITKCommon ITKIO ITKNumerics)

add_executable(procgems procgems.cpp)
target_link_libraries(procgems FID procparser cmdparser ITKCommon ITKIO ITKNumerics)

add_executable(procfsems procfsems.cpp)
target_link_libraries(procfsems FID procparser cmdparser ITKCommon ITKIO ITKNumerics)

add_executable(proc3d proc3d.cpp)
target_link_libraries(proc3d FID procparser cmdparser ITKCommon ITKIO ITKNumerics)

add_executable(procepi procepi.cpp)
target_link_libraries(procepi FID procparser cmdparser ITKCommon ITKIO ITKNumerics)

add_executable(procmsme procmsme.cpp)
target_link_libraries(procmsme FID procparser cmdparser mrfit ITKCommon ITKIO ITKNumerics)

add_executable(procmems procmems.cpp)
target_link_libraries(procmems FID procparser cmdparser mrfit ITKCommon ITKIO ITKNumerics)

add_executable(procmsepi procmsepi.cpp)
target_link_libraries(procmsepi FID procparser cmdparser dtifit ITKCommon ITKIO ITKNumerics)

#add_executable(proclooklocker proclooklocker.cpp)
#target_link_libraries(proclooklocker FID procparser cmdparser mrfit ITKCommon ITKIO ITKNumerics)

