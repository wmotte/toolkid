cmake_minimum_required(VERSION 2.6)
project( resampleDTI )

include_directories( ../cmdparser )

# resample dti...
add_executable( resample_dti 
  ResampleDTI.cxx 
  itkDiffusionTensor3DAbsCorrection.h
  itkDiffusionTensor3DAffineTransform.h
  itkDiffusionTensor3DBSplineInterpolateImageFunction.h
  itkDiffusionTensor3DExtended.h
  itkDiffusionTensor3DFSAffineTransform.h
  itkDiffusionTensor3DInterpolateImageFunction.h
  itkDiffusionTensor3DInterpolateImageFunctionReimplementation.h
  itkDiffusionTensor3DLinearInterpolateFunction.h
  itkDiffusionTensor3DMatrix3x3Transform.h
  itkDiffusionTensor3DNearestCorrection.h
  itkDiffusionTensor3DNearestNeighborInterpolateFunction.h
  itkDiffusionTensor3DNearestNeighborInterpolateFunction.txx
  itkDiffusionTensor3DNonRigidTransform.h
  itkDiffusionTensor3DPPDAffineTransform.h
  itkDiffusionTensor3DRead.h
  itkDiffusionTensor3DResample.h
  itkDiffusionTensor3DRigidTransform.h
  itkDiffusionTensor3DTransform.h
  itkDiffusionTensor3DWindowedSincInterpolateImageFunction.h
  itkDiffusionTensor3DWrite.h
  itkDiffusionTensor3DZeroCorrection.h
  itkMatrixExtended.h
  itkDiffusionTensor3DAffineTransform.txx
  itkDiffusionTensor3DBSplineInterpolateImageFunction.txx
  itkDiffusionTensor3DExtended.txx
  itkDiffusionTensor3DFSAffineTransform.txx
  itkDiffusionTensor3DInterpolateImageFunctionReimplementation.txx
  itkDiffusionTensor3DInterpolateImageFunction.txx
  itkDiffusionTensor3DLinearInterpolateFunction.txx
  itkDiffusionTensor3DMatrix3x3Transform.txx
  itkDiffusionTensor3DNearestNeighborInterpolateFunction.txx
  itkDiffusionTensor3DNonRigidTransform.txx
  itkDiffusionTensor3DPPDAffineTransform.txx
  itkDiffusionTensor3DRead.txx
  itkDiffusionTensor3DResample.txx
  itkDiffusionTensor3DRigidTransform.txx
  itkDiffusionTensor3DTransform.txx
  itkDiffusionTensor3DWindowedSincInterpolateImageFunction.txx
  itkDiffusionTensor3DWrite.txx
  itkMatrixExtended.txx
)

target_link_libraries( resample_dti cmdparser ITKBasicFilters ITKCommon ITKIO )


