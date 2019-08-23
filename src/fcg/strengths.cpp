#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageSpatialObject.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "tkdCmdParser.h"
#include "flens.h"
#include "itkConvertPixelBuffer.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>
#include <cmath>

/**
 Author: W.M. Otte
 Date: 09-06-2009
 Aim: implement strengths_und.m and strengths_dir.m (brain connectivity toolbox Olaf Sporns).
 Strenghts is similar to degree (applicable for binary graphs only).

 Strength is the sum of all connection weights for individual nodes.
 In a directed graph, instrength (outstrength) is the sum of incoming (outgoing) connection weights for individual nodes.
 The strength is the sum of instrength and outstrength.

 Reference: Barrat A, Barthelemy M, Pastor-Satorras R, Vespignani A (2004):
 The architecture of complex weighted networks. Proc Natl Acad Sci U S A 101:3747-3752.
 */

template <typename T>
struct SquareOfDiff:
	public std::binary_function<T, double, T>
{
	T operator()(T x, double mean) const
	{
		x -= (T)mean;
		return x *= x;
	}
};

class Strengths {
public:
	typedef float PixelType;
	typedef unsigned char MaskPixelType;
	typedef itk::Image< MaskPixelType, 3 > MaskImageType;
	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageMaskSpatialObject< 3 > MaskType;
	typedef itk::ImageSpatialObject< 3, PixelType > SpatialType;

	/**
	 * run.
	 */
	void run( const std::string& inputImageName, const std::string& outputFileName, const std::string& maskFileName,
			const std::string& outputImageFileName, float threshold, bool normalize ) {

		// get dimensions...
		unsigned int dims = getDimensions( inputImageName );

		if ( dims == 2 ) {

			// [ 1 ] get strengths and write to txt-output
			std::vector< PixelType > strengths = getStrengths( inputImageName, outputFileName, threshold, normalize );

			if ( ( maskFileName.size() != 0 ) && ( outputImageFileName.size() != 0 ) ) {
				fillImage( maskFileName, outputImageFileName, strengths );
			}

		} else {
			std::cerr << "Number of input (matrix) dimensions should be 2!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

protected:

	/**
	 * Create new image using given mask and fill nonzero voxels with strengths.
	 */
	void fillImage( const std::string& maskFileName, const std::string& outputImageFileName, const std::vector< float >& strengths ) {

		typedef itk::ImageFileReader< MaskImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( maskFileName.c_str() );
		reader -> Update();

		// mask-image (for dims, etc)
		MaskImageType::Pointer maskImage = reader -> GetOutput();
		MaskType::RegionType region = maskImage -> GetLargestPossibleRegion();
		MaskType::SizeType size = region.GetSize();
		unsigned int voxels = region.GetNumberOfPixels();
		const MaskPixelType* dataMask = maskImage->GetPixelContainer()->GetBufferPointer();

		ImageType::Pointer bufferImage = ImageType::New();
		bufferImage->CopyInformation( maskImage );
		bufferImage->SetRegions( region );
		bufferImage->Allocate();
		bufferImage->FillBuffer( 0 );

		float* buffer = bufferImage -> GetPixelContainer() -> GetBufferPointer();

		// -----------------------------------------------------------------------

		// disconnect pipe-lines
		reader = 0;
		//		spatialReader = 0;

		// fill data buffer with strengths
		unsigned int count = 0;
		for ( unsigned int i = 0; i < voxels; ++i ) {
			if ( dataMask[i] != 0 ) {
				buffer[i] = strengths.at( count );
				count++;
			} else {
				buffer[i] = 0;
			}
		}

		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( outputImageFileName.c_str() );
		writer->SetInput( bufferImage );

		try {
			writer->Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << outputImageFileName << std::endl;
		}
	}

	/**
	 * Check dimensions of inputfile. In case of error, EXIT_FAILURE is returned.
	 */
	unsigned int getDimensions( const std::string& inputFileName ) {
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io ) {
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else {
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Get strength process...
	 */
	std::vector< PixelType > getStrengths( const std::string& inputImageName,
									const std::string& outputFileName, float threshold, bool normalize ) {

		typedef itk::Image< PixelType, 2 > ImageType;
		typedef vnl_matrix_ref< PixelType > DataMatrixType;
		typedef float PrecisionType;
		typedef flens::GeMatrix< flens::FullStorage< PrecisionType, flens::RowMajor > > MatrixType;
		typedef vnl_vector< PrecisionType > VectorType;

		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( inputImageName );
		reader -> Update();

		ImageType::Pointer image = reader -> GetOutput();
		reader = 0;

		PixelType* buffer = image -> GetPixelContainer() -> GetBufferPointer();
		ImageType::RegionType region = image -> GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		int rows = size[0];
		int cols = size[1];

		DataMatrixType data( rows, cols, buffer );

		// fill rows with sum of column-values...
		std::vector< PixelType > sums( rows, 0 );

		for ( int i = 0; i < rows; i++ ) {
			sums[i] = 0;
			for ( int j = 0; j < cols; j++ ) {

				if ( data( i,j ) > threshold ) {
					sums[i] += data( i, j );
				}
			}
		}

		// average
		double mean = std::accumulate( sums.begin(), sums.end(), 0.0 ) / sums.size() ;

		// standardize to Z score
		if ( normalize )
		{
			//std::cout << "Standardizing to Z scores" << std::endl;
			// new vector
			std::vector< PixelType > sumsSquared( rows, 0 );
			// fill new vector
			std::copy( sums.begin(), sums.end(), sumsSquared.begin() );
			// square new vector
			std::transform(
				sumsSquared.begin(), sumsSquared.end(),
				sumsSquared.begin(), std::bind2nd( SquareOfDiff<PixelType>(), mean ) );

			PixelType sd = sqrt( std::accumulate( sumsSquared.begin(), sumsSquared.end(), 0.0) / sumsSquared.size() );

			for ( unsigned int i = 0; i < sums.size(); i++ )
			{
				sums[i] = ( sums[i] - mean ) / sd;
			}

			mean = std::accumulate( sums.begin(), sums.end(), 0.0 ) / sums.size() ;
		}

		// write output...
		std::ofstream filestream;
		filestream.open( outputFileName.c_str() );
		if ( filestream.fail() ) {
			std::cerr << "Not able to write to: " << outputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else {
			for ( unsigned int i = 0; i < sums.size(); i++ ) {
				filestream << sums[i] << std::endl;
			}
		}

		// print average strength...
		std::cout << mean << std::endl;

		filestream.close();

		return sums;
	}
};

/**
 * Main: calculate strengths of weighted graphs.
 */
int main( int argc, char ** argv ) {
	tkd::CmdParser p( "strengths", "Output average Strength (sum of sums of all connection weights for individual nodes). \n" );

	std::string inputFileName;
	std::string outputFileName;
	std::string maskFileName;
	std::string outputImageFileName;
	float threshold = 0;
	bool normalize;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Input 2D image: adjacency matrix" ) ->SetRequired(
			true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "output" ) ->SetDescription( "<string> output text file" ) ->SetRequired(
			true );

	p.AddArgument( maskFileName, "mask" ) ->AddAlias( "m" ) ->SetInput( "mask" ) ->SetDescription( "<string> mask image" ) ->SetRequired(
			false );

	p.AddArgument( outputImageFileName, "output-image" ) ->AddAlias( "oi" ) ->SetInput( "output-image" ) ->SetDescription(
			"<string> output image (filled with strengths)" ) ->SetRequired( false );

	p.AddArgument( threshold, "threshold" ) ->AddAlias( "t" ) ->SetInput( "threshold" ) ->SetDescription(
				"<float> threshold (default: 0)" ) ->SetRequired( false );

	p.AddArgument( normalize, "normalize" ) ->AddAlias( "n" ) ->SetInput( "normalize" ) ->SetDescription(
				"<bool> standardize by converting in to Z scores threshold (default: false)" ) ->SetRequired( false );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return -1;
	}

	Strengths strengths;

	strengths.run( inputFileName, outputFileName, maskFileName, outputImageFileName, threshold, normalize );
}

