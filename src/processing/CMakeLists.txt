cmake_minimum_required(VERSION 2.6)
project( processing )


include_directories( ../cmdparser )

# get maximal pixel value coordinates...
add_executable( get_max_dim getMaxDim.cpp )
target_link_libraries( get_max_dim cmdparser ITKCommon ITKIO ITKStatistics )

# diffusion filter...
add_executable( diffusion_filter diffusionFilter.cpp )
target_link_libraries( diffusion_filter cmdparser ITKCommon ITKIO )

# change checkerboard...
add_executable( checkerboard checkerboard.cpp )
target_link_libraries( checkerboard cmdparser ITKCommon ITKIO )

# circularshift
add_executable( circularshift circularshift.cpp )
target_link_libraries( circularshift cmdparser ITKCommon ITKIO )

# convert
add_executable( convert convert.cpp )
target_link_libraries( convert cmdparser ITKCommon ITKIO )

# extract
add_executable( extract extract.cpp )
target_link_libraries( extract cmdparser ITKCommon ITKIO )

# image size
add_executable( image_size imageSize.cpp )
target_link_libraries( image_size cmdparser ITKCommon ITKIO )

# binaryRelabel
add_executable( binary_relabel binaryRelabel.cpp )
target_link_libraries( binary_relabel cmdparser ITKCommon ITKIO ITKNumerics )

# cog4D
add_executable( cog4D cog4D.cpp )
target_link_libraries( cog4D cmdparser ITKCommon ITKIO ITKStatistics )

# cog4D
add_executable( euclidean euclidean.cpp )
target_link_libraries( euclidean cmdparser ITKCommon ITKIO ITKStatistics )

# rescale
add_executable( rescale rescale.cpp )
target_link_libraries( rescale cmdparser ITKCommon ITKIO )

# remove minima
add_executable( remove_minima removeMinima.cpp )
target_link_libraries( remove_minima cmdparser ITKCommon ITKIO )

# mean
add_executable( mean mean.cpp )
target_link_libraries( mean cmdparser ITKCommon ITKIO )

# copy info
add_executable( copyinfo copyinfo.cpp )
target_link_libraries( copyinfo cmdparser ITKCommon ITKIO )

# little2bigendian
add_executable( little2bigendian little2bigendian.cpp itkBigMetaImageIO.h itkBigMetaImageIO.cpp )
target_link_libraries( little2bigendian cmdparser ITKCommon ITKIO )


