#include "tkdCmdParser.h"

#include "itkImage.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "itkImageMaskSpatialObject.h"
#include "vnl/vnl_matrix_ref.h"
#include "nrFourier.h"
#include "graphCommon.h"
#include <vector>
#include <iostream>

/**
 * Segmentation combiner.
 * 
 * For all 4D input label images, construct weighted graph for label similarity.
 *  
 */
class SegmentationCombiner
{

public:

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName, const std::string& maskFileName, const std::string& outputFileName )
	{

		if ( GetDimensions( inputFileName ) == 4 )
		{
			ConstructMatrix( inputFileName, maskFileName, outputFileName );
		} else
		{
			std::cerr << "Input should be 4D!"<< std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	typedef float PixelType;
	typedef float OutputPixelType;
	typedef unsigned char MaskPixelType;

	typedef itk::Image< PixelType, 4 > ImageType;
	typedef itk::Image< MaskPixelType, 3 > MaskImageType;
	typedef itk::ImageMaskSpatialObject< 3 > MaskType;
	typedef vnl_vector< PixelType > VectorType;	
	typedef vnl_matrix_ref< PixelType > DataType;
	typedef vnl_matrix< PixelType > MatrixType;
	typedef itk::Image< PixelType, 3 > MapType;
	typedef vnl_matrix_ref< MaskPixelType > DataMaskType;
	typedef itk::Image< OutputPixelType, 2 > OutputImageType;
	typedef vnl_matrix_ref< OutputPixelType > OutputDataType;

	/**
	 * Load image.
	 */
	ImageType::Pointer loadImage( const std::string& filename )
	{
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New( );
		reader->SetFileName( filename.c_str( ) );
		reader->Update( );

		return reader->GetOutput( );
	}

	/**
	 * Load mask.
	 */
	MaskType::Pointer loadMask( const std::string& filename )
	{
		typedef itk::ImageFileReader< MaskImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New( );
		reader->SetFileName( filename.c_str( ) );
		reader->Update( );

		MaskType::Pointer mask = MaskType::New( );
		mask->SetImage( reader->GetOutput( ) );
		return mask;
	}

	/**
	 * Construct weighted matrix from given 4D input images.
	 */
	void ConstructMatrix( const std::string& inputFileName, 
			const std::string& maskFileName, const std::string& outputFileName )
	{
		std::cout << "Read "<< inputFileName << std::endl;
		ImageType::Pointer image = loadImage( inputFileName );
		std::cout << "Read "<< maskFileName << std::endl;
		MaskType::Pointer mask = loadMask( maskFileName );

		ImageType::RegionType region = image->GetLargestPossibleRegion( );
		ImageType::SizeType size = region.GetSize( );

		int voxels = size[ 0 ] * size[ 1 ]* size[ 2 ];
		int timepoints = size[ 3 ];

		PixelType* data = image->GetPixelContainer()->GetBufferPointer( );
		MaskPixelType* dataMask = const_cast< MaskType::ImageType* >( mask->GetImage() )->GetPixelContainer()->GetBufferPointer( );

		DataType source( timepoints, voxels, data );
		MatrixType matrix = source.extract( source.rows( ), source.cols() );
		
		std::vector< int > indices;
		std::vector< int > inverseIndices;
		int rowCount = 0;
		for ( int i = 0; i < voxels; ++i )
		{
			if ( dataMask[ i ] )
			{
				indices.push_back( i );
				inverseIndices.push_back( rowCount++ );
			} 
			else
			{
				inverseIndices.push_back( -1 );
			}
		}

		int validVoxels = indices.size( );
		std::cout << validVoxels << " voxels with "<< matrix.rows() << " timepoints"<< std::endl;

		OutputImageType::Pointer outputImage = OutputImageType::New( );
		OutputImageType::RegionType outputRegion;
		OutputImageType::SizeType outputSize;
		outputSize[ 0 ] = validVoxels;
		outputSize[ 1 ] = validVoxels;
		outputRegion.SetSize( outputSize );
		outputImage->SetRegions( outputRegion );
		outputImage->Allocate( );
		outputImage->FillBuffer( itk::NumericTraits<OutputPixelType>::Zero );
		OutputPixelType* outputData = outputImage->GetPixelContainer()->GetBufferPointer( );
		OutputDataType output( validVoxels, validVoxels, outputData );

		std::vector< VectorType > vectors;

		for ( int i = 0; i < validVoxels; ++i )
		{
			VectorType a = matrix.get_column( indices[ i ] );
			vectors.push_back( a );
		}

		image = 0;
		
		// first, calculate functional distance ...
		for ( int i = 0; i < validVoxels; ++i )
		{
			const VectorType& a = vectors[ i ];

			for ( int j = i; j < validVoxels; ++j )
			{
				const VectorType& b = vectors[ j ];

				// get raw functional distance ...
				OutputPixelType r = SimilarityWeight( a, b );

				output( i, j ) = r;
				output( j, i ) = r;
			}
		}
	
		// write matrix output...
		typedef itk::ImageFileWriter< OutputImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( outputFileName.c_str() );
		writer->SetInput( outputImage );

		std::cout << "Write " << outputFileName << std::endl;
		writer->Update();
	}

	/**
	 * Return weight of given integer vectors between 0 and 1.
	 */
	PixelType SimilarityWeight( const VectorType& a, const VectorType& b )
	{
		if( a.size() != b.size() )
		{
			std::cerr << "*** ERROR ***: could not calculate similarity weight (vectors not equal in size)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		unsigned int weight = 0;
		for( unsigned int i = 0; i < a.size(); i++ )
		{
			if( a[i] == b[i] )
			weight++;
		}
		return weight / static_cast< PixelType >( a.size() );
	}

	/**
	 * Check dimensions of inputfile. In case of error, EXIT_FAILURE is returned.
	 */
	unsigned int GetDimensions( const std::string& inputFileName )
	{
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io )
		{
			std::cerr << "Could not read: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else
		{
			io->SetFileName( inputFileName );
			io->ReadImageInformation();
			return io->GetNumberOfDimensions();
		}
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "segmentation combiner", "Combine segmentation labels in 4D input into a weighted similarity graph");

	std::string inputFileName;
	std::string maskFileName;
	std::string outputFileName;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Input 4D image: aligned 3D segmentations" ) ->SetRequired( true );

	p.AddArgument( maskFileName, "mask" ) ->AddAlias( "m" ) ->SetInput( "filename" ) ->SetDescription( "Mask for voxels of interest" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "float" ) ->SetDescription( "Weighted output matrix" ) ->SetRequired( true );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	SegmentationCombiner segmentationCombiner;

	segmentationCombiner.Run( inputFileName, maskFileName, outputFileName );

	return EXIT_SUCCESS;
}

