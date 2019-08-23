#ifndef __Relabel_h__
#define __Relabel_h__

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkLabelGeometryImageFilter.h"

static const unsigned int Dimension = 3;
typedef unsigned short LabelPixelType;
typedef itk::Image< LabelPixelType, Dimension > LabelImageType;
typedef itk::ImageFileReader< LabelImageType > LabelReaderType;
typedef itk::ImageFileWriter< LabelImageType > LabelWriterType;
typedef itk::RelabelComponentImageFilter< LabelImageType, LabelImageType > RelabelType;
typedef itk::LabelGeometryImageFilter< LabelImageType > LabelGeometryType;
typedef LabelImageType::Pointer LabelImagePointerType;

typedef double PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef ImageType::Pointer ImagePointerType;

/**
 * Relabel image.
 */
class Relabel
{
public:

	void Run( const std::string& inputFileName, const std::string& refFileName, const std::string& outputFileName );
	void RelabelImage( LabelImagePointerType& output, const LabelImagePointerType& input, const LabelImagePointerType& ref );
	void ReadImage( ImagePointerType& output, const std::string& inputFileName );
	void ConvertImageToLabelImage( LabelImagePointerType& output, ImagePointerType& input );
	void ConvertLabelImageToImage( ImagePointerType& output, LabelImagePointerType& input );
};

#endif /*__Relabel_h__*/
