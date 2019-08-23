#include "annRelabel.h"
#include "tkdCmdParser.h"

/**
 * Run. Relabel 3D input with refNumber the given reference volume.
 */
void Relabel::Run( const std::string& inputFileName, const std::string& refFileName, const std::string& outputFileName )
{
	ImagePointerType inputImage;
	ImagePointerType refImage;
	ReadImage( inputImage, inputFileName );
	ReadImage( refImage, refFileName );

	LabelImagePointerType labelInputImage;
	LabelImagePointerType labelRefImage;

	ConvertImageToLabelImage( labelInputImage, inputImage );
	ConvertImageToLabelImage( labelRefImage, refImage );

	LabelImagePointerType labelOutputImage;
	RelabelImage( labelOutputImage, labelInputImage, labelRefImage );

	// Write labelOutputImage to disk...
	LabelWriterType::Pointer writer = LabelWriterType::New();
	writer->SetFileName( outputFileName.c_str() );
	writer->SetInput( labelOutputImage );
	writer->Update();
}

/**
 * Find closest labels with ref.
 */
void Relabel::RelabelImage( LabelImagePointerType& output, const LabelImagePointerType& input, const LabelImagePointerType& ref )
{
	// image
	LabelGeometryType::Pointer imageGeometryFilter = LabelGeometryType::New();
	imageGeometryFilter->SetInput( input );
	imageGeometryFilter->Update();
	LabelGeometryType::LabelsType imageLabels = imageGeometryFilter->GetLabels();
	LabelGeometryType::LabelsType::iterator imageIt;

	// ref
	LabelGeometryType::Pointer refGeometryFilter = LabelGeometryType::New();
	refGeometryFilter->SetInput( ref );
	refGeometryFilter->Update();
	LabelGeometryType::LabelsType refLabels = refGeometryFilter->GetLabels();
	LabelGeometryType::LabelsType::iterator refIt;

	std::vector< unsigned int > newValues( imageGeometryFilter->GetNumberOfLabels() );

	if( refGeometryFilter->GetNumberOfLabels() != imageGeometryFilter->GetNumberOfLabels() )
	{
		std::cerr << "*** ERROR ***: number of labels in reference is not equal to number of labels in image!" << std::endl;
		exit( EXIT_FAILURE );
	}

	for ( imageIt = imageLabels.begin(); imageIt != imageLabels.end(); ++imageIt )
	{
		LabelImageType::PointType imCentroid = imageGeometryFilter->GetCentroid( *imageIt );
		PixelType minDistance = itk::NumericTraits< PixelType >::max();

		if ( *imageIt != 0 )
		{
			for ( refIt = refLabels.begin(); refIt != refLabels.end(); ++refIt )
			{
				if ( *refIt != 0 )
				{
					LabelImageType::PointType rCentroid = refGeometryFilter->GetCentroid( *refIt );
					PixelType distance = imCentroid.EuclideanDistanceTo( rCentroid );

					if ( distance < minDistance )
					{
						minDistance = distance;
						newValues.at( *imageIt ) = *refIt;
					}
				}
			}
		}
	}

	// Create new, empty image, and iterate over Regions -> remap
	output = LabelImageType::New();
	output->CopyInformation( input );
	output->SetRegions( input->GetLargestPossibleRegion() );
	output->Allocate();
	output->FillBuffer( 0 );

	itk::ImageRegionIterator< LabelImageType > it( input, input->GetLargestPossibleRegion() );
	itk::ImageRegionIterator< LabelImageType > nit( output, output->GetLargestPossibleRegion() );

	while( !it.IsAtEnd() )
	{
		if( !it.Get() == 0 )
		{
			nit.Set( newValues.at( it.Get() ) );
		}
		++it;
		++nit;
	}
}

/**
 * Read 3D double image.
 */
void Relabel::ReadImage( ImagePointerType& output, const std::string& inputFileName )
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inputFileName );
	reader->Update();
	output = reader->GetOutput();
}

/**
 * Convert 3D double image to label 3D image.
 */
void Relabel::ConvertImageToLabelImage( LabelImagePointerType& output, ImagePointerType& input )
{
	typedef itk::CastImageFilter< ImageType, LabelImageType > CastImageFilterType;
	CastImageFilterType::Pointer filter = CastImageFilterType::New();
	filter->SetInput( input );
	filter->Update();
	output = filter->GetOutput();
}

/**
 * Convert 3D label image to double 3D image.
 */
void Relabel::ConvertLabelImageToImage( ImagePointerType& output, LabelImagePointerType& input )
{
	typedef itk::CastImageFilter< LabelImageType, ImageType > CastImageFilterType;
	CastImageFilterType::Pointer filter = CastImageFilterType::New();
	filter->SetInput( input );
	filter->Update();
	output = filter->GetOutput();
}
