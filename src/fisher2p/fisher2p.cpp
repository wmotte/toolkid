#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include "tkdCmdParser.h"

#include <boost/math/special_functions/erf.hpp>


/**
 * Fisher z' to p-value (Numerical Recipes, 746)
 */
class Fisher2P
{
	typedef double PixelType;
	typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageReaderType;
	typedef itk::ImageFileWriter< ImageType > ImageWriterType;
    typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
    typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ConstIteratorType;

	public:

    void Run( const std::string& inputFileName, unsigned int N )
	{
		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName( inputFileName.c_str() );
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();

        ImageType::Pointer negative = ImageType::New();
        negative->SetRegions( image->GetLargestPossibleRegion() );
        negative->SetSpacing( image->GetSpacing() );
        negative->SetOrigin( image->GetOrigin() );
        //negative->CopyInformation( image );
        negative->Allocate();
        negative->FillBuffer( 0 );

        ImageType::Pointer positive = ImageType::New();
        positive->SetRegions( image->GetLargestPossibleRegion() );
        positive->SetSpacing( image->GetSpacing() );
        positive->SetOrigin( image->GetOrigin() );
        positive->Allocate();
        positive->FillBuffer( 0 );

        ConstIteratorType it( image, image->GetLargestPossibleRegion() );
        it.GoToBegin();

        while( ! it.IsAtEnd() )
        {
            ImageType::IndexType index = it.GetIndex();
            PixelType z = image->GetPixel( index );
            PixelType p = P( z, N );

            std::cout << "P-value: " << p << std::endl;

            if( z < 0 )
                negative->SetPixel( index, p );
            else
                positive->SetPixel( index, p );
            
            ++it;
        }

        // Write
        std::string outputFileNameNegative = "neg.nii.gz";
        std::string outputFileNamePositive = "pos.nii.gz";

        if( !outputFileNameNegative.empty() )
        {
            std::cout << "Write p-values for negative z' values to: " << outputFileNameNegative << std::endl;
            ImageWriterType::Pointer writer = ImageWriterType::New();
            writer->SetInput( negative );
            writer->SetFileName( outputFileNameNegative );
            writer->Update();
        }
        if( !outputFileNamePositive.empty() )
        {
            std::cout << "Write p-values for positive z' values to: " << outputFileNamePositive << std::endl;
            ImageWriterType::Pointer writer = ImageWriterType::New();
            writer->SetInput( positive );
            writer->SetFileName( outputFileNamePositive );
            writer->Update();
        }
	}


	private:

    double P( double z, double N )
    {
        if( N < 4 )
        {
            std::cerr << "*** ERROR ***: N should be larger than 3!" << std::endl;
            exit( EXIT_FAILURE );
        }
        // Like randomise output ...
        //return 1 - erfc( std::abs( z ) * std::sqrt( N - 3 ) ) / static_cast< double >( sqrt( 2. ) );
        return 1 - erfc( std::abs( z ) * std::sqrt( N - 3 ) ) / static_cast< double >( sqrt( 2. ) );
    }

};


int main( int argc, char ** argv )
{
    tkd::CmdParser p( "fisher2p", "Calculate p-value from given z' value and N" );
    std::string inputFileName;
    int N;

    p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "input" )
    ->SetDescription( "Fisher z' input file name" )
    ->SetRequired( true );

    p.AddArgument( N, "N" )
    ->SetInput( "N" )
    ->SetDescription( "Sample size used to calculate z'" )
    ->SetRequired( true );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  Fisher2P fisher2p;
  fisher2p.Run( inputFileName, N );

  return EXIT_SUCCESS;
}
