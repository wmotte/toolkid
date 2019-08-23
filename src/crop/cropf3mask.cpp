#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageMaskSpatialObject.h"

int main(int argc, char **argv)
{
	if (argc != 4)
	{
		std::cout << "Usage: " << argv[0] << " input mask output" << std::endl;
		return -1;
	}
	
	typedef float PixelType;
	typedef unsigned char MaskPixelType;
	typedef itk::Image<PixelType, 3> ImageType;
	typedef itk::Image<MaskPixelType, 3> MaskImageType;
	typedef itk::ImageMaskSpatialObject<3> MaskType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;
	typedef itk::ExtractImageFilter<ImageType, ImageType> FilterType;

	ReaderType::Pointer reader = ReaderType::New();
	MaskReaderType::Pointer maskReader = MaskReaderType::New();
	FilterType::Pointer filter = FilterType::New();
	WriterType::Pointer writer = WriterType::New();	
	reader->SetFileName(argv[1]);
	maskReader->SetFileName(argv[2]);
	filter->SetInput(reader->GetOutput());
	writer->SetInput(filter->GetOutput());
	writer->SetFileName(argv[3]);
	
	maskReader->Update();
	MaskImageType::Pointer maskImage = maskReader->GetOutput();
	MaskType::Pointer mask = MaskType::New();
	mask->SetImage(maskImage);
	
	MaskType::RegionType region = mask->GetAxisAlignedBoundingBoxRegion();
	
	filter->SetExtractionRegion(region);
	writer->Update();
	
	return 0;
}
