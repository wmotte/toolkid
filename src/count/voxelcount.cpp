#include "tkdCmdParser.h"

#include "graphCommon.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImage.h"

/**
 * wim@invivonmr.uu.nl, 26-08-2010.
 * VoxelCount voxels above threshold.
 * Sum of voxels (weighted, above threshold, 7-10-2010
 */
namespace count
{
	typedef double PixelType;

	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;


	// *******************************************************************************

	/**
	 * Option list.
	 */
	struct parameters
	{
		std::string inputFileName;
		std::string maskFileName;
        PixelType threshold;
        bool weighted;
        bool average;
	};

	// *******************************************************************************

	/**
	 *
	 */
	class VoxelCount
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
			ImageType::Pointer image = ReadImage( args.inputFileName );

			double count = 0;

			if( args.maskFileName.empty() )
			{
				// Create datacontainer.
				count = CountVoxels( image, args.threshold, args.weighted, args.average );
			}
			else
			{
				ImageType::Pointer mask = ReadImage( args.maskFileName );
				count = CountVoxelsWithinMask( image, mask, args.threshold, args.weighted, args.average );
			}

            
            // Output.
            std::cout << count << std::endl;

			exit( EXIT_SUCCESS );
		}

	private:

		/**
		 * VoxelCount voxels above threshold.
		 */
		double CountVoxelsWithinMask( const ImageType::Pointer& image, const ImageType::Pointer& mask, double threshold, bool weighted, bool average )
		{
		    IteratorType it( image, image->GetLargestPossibleRegion() );
            IteratorType mit( mask, mask->GetLargestPossibleRegion() );
		    double count = 0;
            double totalvoxels = 0;
			for ( mit.GoToBegin(), it.GoToBegin(); ! it.IsAtEnd(), ! mit.IsAtEnd(); ++it, ++mit )
			{
			    PixelType takePixel = image->GetPixel( it.GetIndex() );

			    if( mask->GetPixel( it.GetIndex() ) > 0 )
			    {
			    	if ( takePixel >= threshold )
			    	{
                        totalvoxels++;
			    		if( weighted )
                            count += takePixel;
                        else
                            count += 1;
			    	}
			    }
		    }
            if( average )
                return count / totalvoxels;
            else
                return count;
		}

		/**
		 * VoxelCount voxels above threshold.
		 */
		double CountVoxels( const ImageType::Pointer& image, double threshold, bool weighted, bool average )
		{
		    IteratorType it( image, image->GetLargestPossibleRegion() );
            double count = 0;
            double totalvoxels = 0;
			for ( it.GoToBegin(); ! it.IsAtEnd(); ++it )
			{
			    PixelType takePixel = image->GetPixel( it.GetIndex() );

				if ( takePixel >= threshold )
                {
                    totalvoxels++;
				    if( weighted )
                        count += takePixel;
                    else
                	    count += 1;
                }
		    }
            if( average )
                return count / totalvoxels;
            else
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
		 * Checks.
		 */
		void Checks( parameters& args )
		{
			// image dims
			if( graph::Graph< PixelType >::GetImageDimensions( args.inputFileName ) != 3 )
			{
				std::cerr << "*** ERROR ***: input is not 3D!" << std::endl;
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
	tkd::CmdParser p( argv[0], "Return number of voxels >= threshold (or sum if weighted)." );

	count::parameters args;

    args.threshold = 1.0;
    args.weighted = false;

	p.AddArgument( args.inputFileName, "input" )
		->AddAlias( "i" )
		->SetDescription( "3D Input file" )
		->SetMinMax( 1,1 )
		->SetRequired( true );

	p.AddArgument( args.threshold, "threshold" )
		->AddAlias( "t" )
		->SetDescription( "Threshold (default: 1.0)" )
		->SetMinMax( 1, 1 );

	p.AddArgument( args.maskFileName, "mask" )
		->AddAlias( "m" )
		->SetDescription( "Mask (optional)" )
		->SetMinMax( 1, 1 );
	
    p.AddArgument( args.weighted, "weighted" )
		->AddAlias( "w" )
		->SetDescription( "Sum voxel values in place of count (default: false)" )
		->SetMinMax( 1, 1 );
    
    p.AddArgument( args.average, "average" )
		->AddAlias( "a" )
		->SetDescription( "Average voxel values in place of count (default: false)" )
		->SetMinMax( 1, 1 );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	count::VoxelCount count;

	count.Run( args );

	return EXIT_SUCCESS;
}
