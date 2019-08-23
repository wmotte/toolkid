#ifndef __rbmImage_h__
#define __rbmImage_h__
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "vnl/vnl_matrix_ref.h"
#include "vtkSmartPointer.h"
#include "itkLightObject.h"

class vtkImageData;

namespace rbm
{

class Image : public itk::LightObject
{
public:
	typedef Image Self;
	typedef itk::LightObject Superclass;
	typedef itk::SmartPointer< Self > Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	itkNewMacro( Self );

	typedef float PixelType;
	typedef itk::Image< PixelType, 4 > SeriesType;
	typedef itk::Image< PixelType, 3 > ImageType;
	typedef SeriesType::RegionType RegionType;
	typedef SeriesType::SizeType SizeType;
	typedef SeriesType::IndexType VoxelType;
	typedef SeriesType::PointType CoordinateType;
	typedef SeriesType::SpacingType SpacingType;
	typedef itk::ImageToVTKImageFilter< ImageType > FilterType;
	typedef itk::ExtractImageFilter< SeriesType, ImageType > ExtractType;
	typedef vnl_matrix_ref< PixelType > MatrixType;
	typedef vnl_vector< PixelType > VectorType;

	void SetSeries( SeriesType::Pointer series, int volume = 0 );

	ImageType::Pointer GetImage() const;
	vtkSmartPointer< vtkImageData > GetImageData() const;
	SeriesType::Pointer GetSeries() const;
	MatrixType GetImageMatrix() const;
	SpacingType GetSpacing() const;
	CoordinateType GetOrigin() const;

	int GetNumberOfVolumes() const;
	int GetNumberOfVoxels() const;
	SizeType GetSize() const;
	void SetVolume( int volume );
	int GetVolume() const;
	double GetZ() const;

	VoxelType CoordinateToVoxel( const CoordinateType& coordinate ) const;
	CoordinateType VoxelToCoordinate( const VoxelType& voxel ) const;

	int GetProfile( VoxelType start, VoxelType end, VectorType& mean, VectorType& sd ) const;

protected:
	Image();

	void Initialize( SeriesType::Pointer series, int volume = 0 );
	void Update();

	SeriesType::Pointer m_Series;
	ImageType::Pointer m_Image;
	FilterType::Pointer m_Filter;
	ExtractType::Pointer m_Extract;
	vtkSmartPointer< vtkImageData > m_ImageData;
	int m_Volume;
	double m_Z;
	int m_Volumes;
	int m_Voxels;
	RegionType m_Region;
	SizeType m_Size;
	RegionType m_ExtractionRegion;
	SizeType m_ExtractionSize;
	VoxelType m_ExtractionIndex;
	CoordinateType m_Origin;
	SpacingType m_Spacing;
	MatrixType m_Matrix;

private:
	Image( const Image& other ); // not implemented
	void operator=( const Image& other ); // not implemented
};

} // end namespace rbm

#endif /*__rbmImage_h__*/
