
#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "Crop.h"
/**
 * In this class the input image will be cropped automatically, between all voxels
 * with value < minValue.
 */


	/**
	 * Process.
	 */
	template< int dimension >
	void CropTemplate<dimension>::run ( const std::string& input, const std::string& output, float minValue ) {

		//const unsigned int Dimension = 3;

		typedef itk::Image< float, Dimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;
		typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractImageFilterType;

		// get input image...
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		ImageType::ConstPointer image = reader -> GetOutput();
		reader -> Update();

		// find index and size to crop...
		IteratorType it( image, image -> GetLargestPossibleRegion() );

		ImageType::IndexType startCropIndex = image -> GetLargestPossibleRegion().GetIndex();
		ImageType::IndexType endCropIndex = image -> GetLargestPossibleRegion().GetIndex();
		ImageType::SizeType outputSize = image -> GetLargestPossibleRegion().GetSize();

		// set start crop index to maximum...
		for ( unsigned int i = 0; i < Dimension; i++ ) {
			startCropIndex[i] = outputSize[i] - 1;
		}

		// loop over voxels. If voxel value > minimum value, update start and end crop indices...
		for ( ; !it.IsAtEnd(); ++it ) {

			ImageType::IndexType index = it.GetIndex();
			float intensity = image -> GetPixel( index );

			if ( intensity > minValue ) {
				for ( unsigned int i = 0; i < Dimension; i++ ) {
					if ( endCropIndex[i] < index[i] ) {
						endCropIndex[i] = index[i];
					} else if ( startCropIndex[i] > index[i] ) {
						startCropIndex[i] = index[i];
					}
				}
			}
		}

		// set output size...
		for ( unsigned int i = 0; i < Dimension; i++ ) {
			outputSize[i] = ( endCropIndex[i] - startCropIndex[i] ) + 1;
		}

		// set crop region...
		ImageType::RegionType outputRegion;
		outputRegion.SetSize( outputSize );
		outputRegion.SetIndex( startCropIndex );

		// filter...
		ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
		extractFilter -> SetInput( image );
		extractFilter -> SetExtractionRegion( outputRegion );

		// write output...
		WriterType::Pointer writer = WriterType::New();
		writer -> SetFileName( output );
		writer -> SetInput( extractFilter -> GetOutput() );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}
	}
};
