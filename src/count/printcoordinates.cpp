#include "tkdCmdParser.h"

#include "graphCommon.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

/**
 * 14-09-2010.
 * Get voxel from central mask and compute for all voxels within mask the relative x,y,z coordinates
 */
namespace count
{
	typedef double PixelType;

	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
	typedef itk::ImageFileWriter< ImageType > ImageFileWriterType;
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;


	// *******************************************************************************

	/**
	 * Option list.
	 */
	struct parameters
	{
		std::string centerFileName;
		std::string maskFileName;
		std::string outputXFileName;
		std::string outputYFileName;
		std::string outputZFileName;

	};

	// *******************************************************************************

	/**
	 *
	 */
	class PrintCoordinates
	{

	public:

		/**
		 * Run.
		 */
		void Run( parameters& args )
		{
			// Check args.
			Checks( args );

			// Read input.
			ImageType::Pointer mask = ReadImage( args.maskFileName );
			ImageType::Pointer center = ReadImage( args.centerFileName );

			ImageType::IndexType index = GetCenterVoxelIndex( center );

			ImageType::Pointer outputX = ImageType::New();
			ImageType::Pointer outputY = ImageType::New();
			ImageType::Pointer outputZ = ImageType::New();

			GetCoordinates( index, mask, outputX, outputY, outputZ );

			WriteImage( args.outputXFileName, outputX );
			WriteImage( args.outputYFileName, outputY );
			WriteImage( args.outputZFileName, outputZ );

			exit( EXIT_SUCCESS );
		}

	private:

		/**
		 * Get coordinates relative to index voxel.
		 */
		void GetCoordinates( const ImageType::IndexType& index, const ImageType::Pointer& mask,
				ImageType::Pointer& outputX,ImageType::Pointer& outputY, ImageType::Pointer& outputZ )
		{
			outputX->CopyInformation( mask );
			outputX->SetRegions( mask->GetLargestPossibleRegion() );
			outputX->Allocate();
			outputX->FillBuffer( 0 );

			outputY->CopyInformation( mask );
			outputY->SetRegions( mask->GetLargestPossibleRegion() );
			outputY->Allocate();
			outputY->FillBuffer( 0 );

			outputZ->CopyInformation( mask );
			outputZ->SetRegions( mask->GetLargestPossibleRegion() );
			outputZ->Allocate();
			outputZ->FillBuffer( 0 );

			IteratorType it( mask, mask->GetLargestPossibleRegion() );

			for ( it.GoToBegin(); ! it.IsAtEnd(); ++it )
			{
				if( it.Get() > 0 )
				{
					ImageType::IndexType pixelIndex = it.GetIndex();

					outputX->SetPixel( pixelIndex, pixelIndex[0] - index[0] );
					outputY->SetPixel( pixelIndex, pixelIndex[1] - index[1] );
					outputZ->SetPixel( pixelIndex, pixelIndex[2] - index[2] );

				}
			}
		}

		/**
		 * Return point of first non-zero voxel.
		 */
		ImageType::IndexType GetCenterVoxelIndex( const ImageType::Pointer& mask )
		{
			IteratorType mit( mask, mask->GetLargestPossibleRegion() );
			for ( mit.GoToBegin(); ! mit.IsAtEnd(); ++mit )
				if( mask->GetPixel( mit.GetIndex() ) > 0 )
					return mit.GetIndex();

			std::cerr << "*** ERROR ***: no center voxel found!" << std::endl;
			exit( EXIT_FAILURE );
		}

		/**
		 * VoxelCount voxels above threshold.
		 */
		unsigned int CountVoxels( const ImageType::Pointer& image, double threshold )
		{
		    IteratorType it( image, image->GetLargestPossibleRegion() );
            unsigned int count = 0;

			for ( it.GoToBegin(); ! it.IsAtEnd(); ++it )
			{
			    PixelType takePixel = image->GetPixel( it.GetIndex() );

				if ( takePixel >= threshold )
							count++;
		    }
            return count;
		}

		/**
		 * Read image.
		 */
		ImageType::Pointer ReadImage( const std::string& image )
		{
			ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
			reader->SetFileName( image );
			reader->Update();
			return reader->GetOutput();
		}

		/**
		 * Write.
		 */
		void WriteImage( const std::string& outputFileName, const ImageType::Pointer& image )
		{
			ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
			writer->SetFileName( outputFileName );
			writer->SetInput( image );

			try
			{
				writer->Update();
			}
			catch( itk::ExceptionObject& e )
			{
				std::cerr << "*** ERROR ***, could not write to: " << outputFileName << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		/**
		 * Checks.
		 */
		void Checks( parameters& args )
		{
			// image dims
			if( graph::Graph< PixelType >::GetImageDimensions( args.maskFileName ) != 3 )
			{
				std::cerr << "*** ERROR ***: mask is not 3D!" << std::endl;
				exit( EXIT_FAILURE );
			}
			if( graph::Graph< PixelType >::GetImageDimensions( args.centerFileName ) != 3 )
			{
				std::cerr << "*** ERROR ***: center is not 3D!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}
	};

} // end namespace count

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Calculate x,y,z coordinates relative to center voxel." );

	count::parameters args;

	p.AddArgument( args.maskFileName, "mask" )
		->AddAlias( "m" )
		->SetDescription( "3D Mask file" )
		->SetMinMax( 1,1 )
		->SetRequired( true );

	p.AddArgument( args.centerFileName, "center" )
		->AddAlias( "c" )
		->SetDescription( "3D Mask file with one voxel taken as center coordinate" )
		->SetMinMax( 1,1 )
		->SetRequired( true );

	p.AddArgument( args.outputXFileName, "output-x" )
		->AddAlias( "ox" )
		->SetDescription( "3D output file, relative x coordinate" )
		->SetMinMax( 1,1 )
		->SetRequired( true );

	p.AddArgument( args.outputYFileName, "output-y" )
		->AddAlias( "oy" )
		->SetDescription( "3D output file, relative y coordinate" )
		->SetMinMax( 1,1 )
		->SetRequired( true );

	p.AddArgument( args.outputZFileName, "output-z" )
		->AddAlias( "oz" )
		->SetDescription( "3D output file, relative z coordinate" )
		->SetMinMax( 1,1 )
		->SetRequired( true );


	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	count::PrintCoordinates app;

	app.Run( args );

	return EXIT_SUCCESS;
}
