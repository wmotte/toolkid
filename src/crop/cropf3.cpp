#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"

int main(int argc, char **argv)
{
	if (argc != 9)
	{
		std::cout << "Usage: " << argv[0] << " input output x y z width height depth" << std::endl;
		return -1;
	}
	
	typedef float PixelType;
	typedef itk::Image<PixelType, 3> ImageType;

	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;
	typedef itk::ExtractImageFilter<ImageType, ImageType> FilterType;

	ReaderType::Pointer reader = ReaderType::New();
	FilterType::Pointer filter = FilterType::New();
	WriterType::Pointer writer = WriterType::New();	
	reader->SetFileName(argv[1]);
	filter->SetInput(reader->GetOutput());
	writer->SetInput(filter->GetOutput());
	writer->SetFileName(argv[2]);
	
	ImageType::RegionType region;
	ImageType::SizeType size;
	ImageType::IndexType index;

	index[0] = atoi(argv[3]);
	index[1] = atoi(argv[4]);
	index[2] = atoi(argv[5]);

	size[0] = atoi(argv[6]);
	size[1] = atoi(argv[7]);
	size[2] = atoi(argv[8]);

	region.SetSize(size);
	region.SetIndex(index);
		
	filter->SetExtractionRegion(region);
	writer->Update();
	
	return 0;
}
