#include "tkdCmdParser.h"

#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "itkImage.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"

#include "graphCommon.h"
#include "annCommon.h"
#include "annRelabel.h"

/**
 * Circular bootstrap.
 *
 * Neuroimage (2010) Pierre Bellec et al,
 * "Multi-level bootstrap analysis of stable cluster in resting-state fMRI".
 *
 * wim@invivonmr.uu.nl, 06-04-2010.
 *
 */
namespace bootstrap
{
	typedef float PixelType;

	typedef std::vector< unsigned int > VectorType;

	typedef itk::Image< PixelType, 3 > Image3DType;
	typedef itk::Image< PixelType, 4 > Image4DType;

	typedef itk::ImageRegionConstIteratorWithIndex< Image4DType > ConstIterator4DType;
	typedef itk::ImageRegionIteratorWithIndex< Image4DType > Iterator4DType;
	typedef itk::ExtractImageFilter< Image4DType, Image4DType > ExtractImageFilterType;
	typedef itk::ImageFileWriter< Image4DType > ImageFileWriterType;

	typedef boost::mt19937 random_number_type;
	typedef boost::uniform_int< int > int_distribution_type;
	typedef boost::variate_generator< random_number_type&, int_distribution_type > int_generator_type;

	// *******************************************************************************

	/**
	 * Option list.
	 */
	struct parameters
	{
		std::string inputFileName;
		std::string maskFileName;
		std::string outputFileName;
		int iterations;
    int clusters;
    int repeat;
  };

  // *******************************************************************************

  /**
   * Construct similarity matrix using circular block boostrapping.
   */
  class Bootstrap
  {

    private:

      Image4DType::Pointer m_Data;
      Image4DType::Pointer m_DoubleData;
      Image3DType::Pointer m_Mask;
      std::vector< Image4DType::Pointer > m_BootstrapContainer;

    public:

      /**
       * Run.
       */
      void Run( const parameters& list )
      {
        // [ 1 ] Read data ...
        Read4D( list.inputFileName );
        ReadMask( list.maskFileName );

        // [ 2 ] Get list of random blocksizes and random offsets ...
        random_number_type generator( time( 0 ) );

        // [ 3 ] Construct circular bootstrap 4D images ...
        FillBootstrapContainer( list.iterations, generator );

        // [ 4 ] DEBUG TODO Write all different 4D images to disk...
        typedef itk::ImageFileWriter< Image4DType > ImageFileWriterType;
        ImageFileWriterType::Pointer writer = ImageFileWriterType::New();

        for( unsigned int i = 0; i < m_BootstrapContainer.size(); i++ )
        {
          std::stringstream ss;
          ss << "test_" << i << ".nii.gz";

          writer->SetFileName( ss.str() );
          writer->SetInput( m_BootstrapContainer.at( i ) );
          writer->Update();
        }
      }

    private:

      /**
       * Fill bootstrap container.
       */
      void FillBootstrapContainer( unsigned int n, random_number_type& generator )
      {
        Image4DType::SizeType size = m_Data->GetLargestPossibleRegion().GetSize();
        unsigned int max = size[3];

        for ( unsigned int i = 0; i < n; i++ )
        {
          std::cout << "Bootstrap nr: " << i << std::endl;
          m_BootstrapContainer.push_back( GetBootstrapImage( generator, max ) );
        }
      }

      /**
       * Return bootstrapped 4D image.
       */
      Image4DType::Pointer GetBootstrapImage( random_number_type& generator, unsigned int max )
      {
        // get number of blocks, sufficient to occupy full 4D image ...
        VectorType blockSizes;
        VectorType offsets;

        while ( std::accumulate( blockSizes.begin(), blockSizes.end(), static_cast< unsigned int > ( 0 ) ) < max )
        {
          offsets.push_back( GetRandomNumber( generator, max ) );
          blockSizes.push_back( GetRandomNumber( generator, 3 * sqrt( max ) ) ); // TODO manual set?
        }

        std::vector< Image4DType::Pointer > blocks = GetBlocks( blockSizes, offsets );

        return MergeBlocks( blocks );
      }

      /**
       * Return circular blocks, taken from input images.
       */
      std::vector< Image4DType::Pointer > GetBlocks( const VectorType& blockSizes, const VectorType& offsets )
      {
        std::vector< Image4DType::Pointer > blocks;

        for ( unsigned int i = 0; i < blockSizes.size(); i++ )
        {
          unsigned int blockSize = blockSizes.at( i );
          unsigned int offset = offsets.at( i );

          ExtractImageFilterType::Pointer filter = ExtractImageFilterType::New();
          filter->SetInput( m_DoubleData );

          Image4DType::RegionType region = m_DoubleData->GetLargestPossibleRegion();
          Image4DType::SizeType size = region.GetSize();
          Image4DType::IndexType index = region.GetIndex();

          size[3] = blockSize;
          index[3] = offset;

          region.SetSize( size );
          region.SetIndex( index );

          filter->SetExtractionRegion( region );

          filter->Update();
          blocks.push_back( filter->GetOutput() );

          m_DoubleData->DisconnectPipeline();
        }
        return blocks;
      }

      /**
       * Merge blocks until size is equal to input 4D Image.
       */
      Image4DType::Pointer MergeBlocks( const std::vector< Image4DType::Pointer >& blocks )
      {
        Image4DType::Pointer image = Image4DType::New();

        image->SetRegions( m_DoubleData->GetLargestPossibleRegion() );
        image->SetOrigin( m_DoubleData->GetOrigin() );
        image->SetSpacing( m_DoubleData->GetSpacing() );

        image->Allocate();
        image->FillBuffer( 0 );

        typedef itk::ImageRegionConstIteratorWithIndex< Image4DType > ConstIterator4D;
        typedef itk::ImageRegionIteratorWithIndex< Image4DType > Iterator4D;

        Image4DType::IndexType index = m_DoubleData->GetLargestPossibleRegion().GetIndex();

        for ( unsigned int i = 0; i < blocks.size(); i++ )
        {
          Image4DType::RegionType region = blocks.at( i )->GetLargestPossibleRegion();
          Image4DType::SizeType size = region.GetSize();
          region.SetIndex( index );

          Iterator4D it4( image, region );
          it4.GoToBegin();

          ConstIterator4D it3( blocks.at( i ), blocks.at( i )->GetLargestPossibleRegion() );
          it3.GoToBegin();

          while ( !it3.IsAtEnd(), !it4.IsAtEnd() )
          {
            it4.Set( it3.Get() );
            ++it3;
            ++it4;
          }
          index[ 3 ] += size[ 3 ];
        }

        // Extract first block from image and return ...
        ExtractImageFilterType::Pointer filter = ExtractImageFilterType::New();
        filter->SetInput( image );

        Image4DType::RegionType cropRegion = image->GetLargestPossibleRegion();
        Image4DType::SizeType cropSize = cropRegion.GetSize();
        Image4DType::IndexType cropIndex = cropRegion.GetIndex();

        cropSize[3] = ( m_Data->GetLargestPossibleRegion().GetSize() )[ 3 ];
        cropIndex[3] = 0;

        cropRegion.SetSize( cropSize );
        cropRegion.SetIndex( cropIndex );

        filter->SetExtractionRegion( cropRegion );
        filter->Update();

        return filter->GetOutput();
      }

      /**
       * Return random number.
       */
      int GetRandomNumber( random_number_type& generator, unsigned int max )
      {
        // uniform distribution: 0 -> 4D image size ...
        int_distribution_type int_uni_dist( 0, max - 1 );
        int_generator_type int_distribution( generator, int_uni_dist );

        int value = int_distribution();
        while ( value < 2 )
          value = int_distribution();

        return value;
      }

      /**
       * Read 3D volume.
       */
      void ReadMask( const std::string& inputFileName )
      {
        if ( !inputFileName.empty() )
        {
          typedef itk::ImageFileReader< Image3DType > ReaderType;
          ReaderType::Pointer reader = ReaderType::New();
          reader->SetFileName( inputFileName );
          reader->Update();
          m_Mask = reader->GetOutput();
        } else
        {
          // use 4D input data as reference for full mask...
          Image4DType::SizeType size4D = m_Data->GetLargestPossibleRegion().GetSize();
          Image4DType::IndexType index4D = m_Data->GetLargestPossibleRegion().GetIndex();
          Image3DType::SizeType size3D;
          Image3DType::IndexType index3D;
          size3D[0] = size4D[0];
          size3D[1] = size4D[1];
          size3D[2] = size4D[2];
          index3D[0] = index4D[0];
          index3D[1] = index4D[1];
          index3D[2] = index4D[2];
          Image3DType::RegionType region3D;
          region3D.SetSize( size3D );
          region3D.SetIndex( index3D );
          m_Mask = Image3DType::New();
          m_Mask->SetRegions( region3D );
          m_Mask->Allocate();
          m_Mask->FillBuffer( 1 );
        }
      }

      /**
       * Read 4D data, and set double data object (for wrapping).
       */
      void Read4D( const std::string& inputFileName )
      {
        typedef itk::ImageFileReader< Image4DType > ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( inputFileName );
        reader->Update();
        m_Data = reader->GetOutput();

        Image4DType::RegionType region = m_Data->GetLargestPossibleRegion();
        Image4DType::SizeType size = region.GetSize();
        size[3] *= 2;
        region.SetSize( size );

        m_DoubleData = Image4DType::New();
        m_DoubleData->SetRegions( region );
        m_DoubleData->SetSpacing( m_Data->GetSpacing() );
        m_DoubleData->SetOrigin( m_Data->GetOrigin() );
        m_DoubleData->Allocate();
        m_DoubleData->FillBuffer( 0 );

        ConstIterator4DType it( m_Data, m_Data->GetLargestPossibleRegion() );
        Iterator4DType nit( m_DoubleData, m_DoubleData->GetLargestPossibleRegion() );

        it.GoToBegin();
        nit.GoToBegin();

        // first round ...
        while( !nit.IsAtEnd(), !it.IsAtEnd() )
        {
          nit.Set( it.Get() );
          ++nit;
          ++it;
        }

        // second round ...
        it.GoToBegin();
        while( !nit.IsAtEnd(), !it.IsAtEnd() )
        {
          nit.Set( it.Get() );
          ++nit;
          ++it;
        }
      }
  };
} // end namespace bootstrap

/**
 * Main.
 */
int main( int argc, char ** argv )
{
  tkd::CmdParser p( argv[0], "Circular block boostrap." );

  bootstrap::parameters args;
  args.iterations = 10;

  p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" )->SetDescription( "4D input image" )->SetMinMax( 1, 1 )->SetRequired(
      true );

  p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" )->SetDescription( "Similarity matrix output file" )->SetMinMax( 1, 1 )->SetRequired(
      true );

  p.AddArgument( args.maskFileName, "mask" ) ->AddAlias( "m" )->SetDescription( "Mask 3D image" )->SetMinMax( 1, 1 );

  p.AddArgument( args.iterations, "iterations" )->AddAlias( "n" )->SetDescription( "Bootstrap iterations (default: 10)" )->SetMinMax( 1,
      1 );

  if ( !p.Parse( argc, argv ) )
  {
    p.PrintUsage( std::cout );
    return EXIT_FAILURE;
  }

  bootstrap::Bootstrap boostrapObj;

  boostrapObj.Run( args );

  return EXIT_SUCCESS;
}

