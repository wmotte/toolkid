#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkExtractImageFilter.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_vector.h"
#include "tkdCmdParser.h"

/**
 * Calculate percent change in Cerebral Blood Volume (rCBV)
 *
 * See Mandeville et al., MRM 39 (1998): p. 622
 * Equation A3:
 *
 * V(t) is CBV at time t
 * S(t) is signal at time t
 * Spre is mean signal before contrast agent injection
 * S(0) is mean signal with contrast agent
 * V(0) is base CBV
 *
 * deltaV(t) / V(0) = ln( S(t) / S(0) ) / ln( S(0) / Spre ) )
 *                  = deltaR2star(t) / deltaR2star(0)
 *
 * gives relative CBV.
 */
class RCBV
{
public:
  typedef double PixelType;
  typedef itk::Image< PixelType, 4 > ImageType;
  typedef itk::Image< PixelType, 3 > VolumeType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::ImageFileWriter< VolumeType > VolumeWriterType;
  typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractType;

  RCBV() : m_OutputIndexMin( -1 ), m_OutputIndexMax( -1 )
  {
  }

  void SetOutputRange( int min, int max )
  {
    m_OutputIndexMin = min;
    m_OutputIndexMax = max;
  }

  void Read( const std::vector< std::string >& filenames )
  {
    std::vector< ImageType::Pointer > images;
    int numberOfTimePoints = 0;
    for( std::vector< std::string >::const_iterator i = filenames.begin(); i != filenames.end(); ++i )
      {
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( i->c_str() );
      reader->Update();
      ImageType::Pointer image = reader->GetOutput();
      reader = 0;

      ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
      numberOfTimePoints += size[ 3 ];
      images.push_back( image );
      }

    ImageType::RegionType region = images[ 0 ]->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    size[ 3 ] = numberOfTimePoints;
    region.SetSize( size );

    m_Image = ImageType::New();
    m_Image->CopyInformation( images[ 0 ] );
    m_Image->SetRegions( region );
    m_Image->Allocate();

    ImageType::RegionType outputRegion, inputRegion;
    ImageType::SizeType outputSize, inputSize;
    ImageType::IndexType outputIndex;

    outputIndex.Fill( 0 );

    for( unsigned int i = 0; i < images.size(); ++i )
      {
      inputRegion = images[ i ]->GetLargestPossibleRegion();
      inputSize = inputRegion.GetSize();

      outputSize = inputSize;

      outputRegion.SetSize( outputSize );
      outputRegion.SetIndex( outputIndex );

      itk::ImageRegionIterator< ImageType > itOut( m_Image, outputRegion );
      itk::ImageRegionConstIterator< ImageType > itIn( images[ i ], inputRegion );

      for( ; !itOut.IsAtEnd(); ++itOut, ++itIn )
        {
        itOut.Set( itIn.Value() );
        }

      outputIndex[ 3 ] += inputSize[ 3 ];
      }
  }

  void Write( const std::string& filename )
  {
    ImageType::RegionType region = m_Image->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    ImageType::IndexType index = region.GetIndex();

    if ( m_OutputIndexMin != -1 )
      {
      index[ 3 ] = m_OutputIndexMin;
      size[ 3 ] = m_OutputIndexMax - m_OutputIndexMin;

      std::cout << "Extract " << m_OutputIndexMin << ", " << m_OutputIndexMax << std::endl;
      }

    region.SetSize( size );
    region.SetIndex( index );

    ExtractType::Pointer extract = ExtractType::New();
    extract->SetInput( m_Image );
    extract->SetExtractionRegion( region );

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( extract->GetOutput() );
    writer->Update();
  }

  void WritePre( const std::string& filename )
  {
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( m_Spre );
    writer->Update();
  }

  void WriteS0( const std::string& filename )
  {
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( m_S0 );
    writer->Update();
  }

  void CalculateRCBV()
  {
    ImageType::SizeType size = m_Image->GetLargestPossibleRegion().GetSize();
    int numberOfRows = size[ 3 ];
    vnl_matrix_ref< PixelType > m = vnl_matrix_ref< PixelType >( numberOfRows, m_NumberOfVoxels, m_Image->GetPixelContainer()->GetBufferPointer() );
    vnl_vector_ref< PixelType > Spre = vnl_vector_ref< PixelType >( m_NumberOfVoxels, m_Spre->GetPixelContainer()->GetBufferPointer() );
    vnl_vector_ref< PixelType > S0 = vnl_vector_ref< PixelType >( m_NumberOfVoxels, m_S0->GetPixelContainer()->GetBufferPointer() );

    for( int i = 0; i < m_NumberOfVoxels; ++i )
      {
      for( int t = 0; t < numberOfRows; ++t )
        {
        // rCBV, equation A3 at page 622
        m( t, i ) = vcl_log( m( t, i ) / S0( i ) ) / vcl_log( S0( i ) / Spre( i ) );

        // scale [-1,1] to [100,300]
        m( t, i ) = m( t, i ) * 100 + 100;

        // threshold at -100% and above
        if ( m( t, i ) < 0 )
          {
          m( t, i ) = 0;
          }
        }
      }
  }

  void CalculatePre( int startPre, int endPre, int start0, int end0 )
  {
    ImageType::SizeType size = m_Image->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin = m_Image->GetOrigin();
    ImageType::SpacingType spacing = m_Image->GetSpacing();

    int numberOfRows = size[ 3 ];
    m_NumberOfVoxels = size[ 0 ] * size[ 1 ] * size[ 2 ];

    m_Spre = VolumeType::New();
    VolumeType::RegionType volumeRegion;
    VolumeType::SizeType volumeSize;
    VolumeType::SpacingType volumeSpacing;
    VolumeType::PointType volumeOrigin;

    for( int i = 0; i < 3; ++i )
      {
      volumeSize[ i ] = size[ i ];
      volumeSpacing[ i ] = spacing[ i ];
      volumeOrigin[ i ] = origin[ i ];
      }

    volumeRegion.SetSize( volumeSize );
    m_Spre->SetSpacing( volumeSpacing );
    m_Spre->SetOrigin( volumeOrigin );
    m_Spre->SetRegions( volumeRegion );
    m_Spre->Allocate();

    m_S0 = VolumeType::New();
    m_S0->CopyInformation( m_Spre );
    m_S0->SetRegions( volumeRegion );
    m_S0->Allocate();

    vnl_matrix_ref< PixelType > m = vnl_matrix_ref< PixelType >( numberOfRows, m_NumberOfVoxels, m_Image->GetPixelContainer()->GetBufferPointer() );
    PixelType* outputPre = m_Spre->GetPixelContainer()->GetBufferPointer();

    int numberOfPre = endPre - startPre;
    int numberOf0 = end0 - start0;

    std::cout << "Voxels: " << m_NumberOfVoxels << std::endl;
    std::cout << "Number of time points: " << numberOfRows << std::endl;
    std::cout << "Spre: " << numberOfPre << " [" << startPre << ", " << endPre << ")" << std::endl;
    std::cout << "S0: " << numberOf0 << " [" << start0 << ", " << end0 << ")" << std::endl;

    // calculate mean pre-contrast
    for( int i = 0; i < m_NumberOfVoxels; ++i )
      {
      vnl_vector< PixelType > column = m.get_column( i );

      outputPre[ i ] = 0;

      for( int j = startPre; j < endPre; ++j )
        {
        outputPre[ i ] += column( j );
        }

      outputPre[ i ] /= static_cast< PixelType >( numberOfPre );
      }

    PixelType* output0 = m_S0->GetPixelContainer()->GetBufferPointer();

    // calculate mean post-contrast (from first scans of time series!)
    for( int i = 0; i < m_NumberOfVoxels; ++i )
      {
      vnl_vector< PixelType > column = m.get_column( i );

      output0[ i ] = 0;

      for( int j = start0; j < end0; ++j )
        {
        output0[ i ] += column( j );
        }

      output0[ i ] /= static_cast< PixelType >( numberOf0 );
      }
  }

protected:
  int m_NumberOfVoxels;
  ImageType::Pointer m_Image;
  VolumeType::Pointer m_Spre;
  VolumeType::Pointer m_S0;
  int m_OutputIndexMin;
  int m_OutputIndexMax;
};


int main( int argc, char ** argv )
{

  tkd::CmdParser p( "rcbv", "Calculate rCBV time series image" );

  std::vector< std::string > inputImages;
  std::vector< int > rangeS0, rangePre, rangeOutput;
  std::string outputImage;
  std::string outputPre;
  std::string outputS0;

  p.AddArgument( inputImages, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filenames" )
    ->SetDescription( "Input time-series images" )
    ->SetRequired( true );

  p.AddArgument( rangePre, "range-Spre" )
    ->AddAlias( "rp" )
    ->SetInput( "start end" )
    ->SetDescription( "Index of first and last image in Spre-series" )
    ->SetMinMax( 2, 2)
    ->SetRequired( true );

  p.AddArgument( rangeS0, "range-S0" )
    ->AddAlias( "r" )
    ->SetInput( "start end" )
    ->SetDescription( "Index of first and last image in S0-series" )
    ->SetMinMax( 2, 2 )
    ->SetRequired( true );

  p.AddArgument( rangeOutput, "range-output" )
    ->AddAlias( "ro" )
    ->SetInput( "[start, end)" )
    ->SetDescription( "Index of first and last volume of output image" )
    ->SetMinMax( 2, 2 );

  p.AddArgument( outputImage, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output image series" )
    ->SetRequired( true );

  p.AddArgument( outputPre, "output-Spre" )
    ->AddAlias( "op" )
    ->SetInput( "filename" )
    ->SetDescription( "Output Spre-image series" );

  p.AddArgument( outputS0, "output-S0" )
    ->AddAlias( "o0" )
    ->SetInput( "filename" )
    ->SetDescription( "Output S0-image series" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  RCBV rcbv;

  if ( rangeOutput.size() != 0 )
    {
    rcbv.SetOutputRange( rangeOutput[ 0 ], rangeOutput[ 1 ] );
    }

  rcbv.Read( inputImages );
  rcbv.CalculatePre( rangePre[ 0 ], rangePre[ 1 ], rangeS0[ 0 ], rangeS0[ 1 ] );
  rcbv.CalculateRCBV();
  rcbv.Write( outputImage );

  if ( outputPre != "" )
    {
    rcbv.WritePre( outputPre );
    }

  if ( outputS0 != "" )
    {
    rcbv.WriteS0( outputS0 );
    }

  return 0;
}
