#include "itkOrientImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"
#include <string>
#include "tkdCmdParser.h"

template< class TPixel, unsigned int VDimension >
int OrientImage(
  const std::string& inputFileName,
  const std::string& outputFileName,
  itk::SpatialOrientation::ValidCoordinateOrientationFlags desired = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI,
  itk::SpatialOrientation::ValidCoordinateOrientationFlags given = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID )
{
  typedef itk::Image< TPixel, VDimension > ImageType;
  typedef itk::OrientImageFilter< ImageType, ImageType > FilterType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  typename ReaderType::Pointer reader = ReaderType::New();
  typename FilterType::Pointer filter = FilterType::New();
  typename WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( inputFileName.c_str() );
  writer->SetFileName( outputFileName.c_str() );
  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  reader->Update();

  if ( given == itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID )
    {
    filter->UseImageDirectionOn();
    }
  else
    {
    filter->UseImageDirectionOff();
    filter->SetGivenCoordinateOrientation( given );
    }

  filter->SetDesiredCoordinateOrientation( desired );
  filter->Print( std::cout );

  writer->Update();

  return 0;
}

void ParseOrientation( const std::string& text, itk::SpatialOrientation::ValidCoordinateOrientationFlags& orient )
{
#define poMacro( t ) \
  if ( text == #t ) \
    { \
    orient = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_##t; \
    return; \
    }

  poMacro( RIP );
  poMacro( LIP );
  poMacro( RSP );
  poMacro( LSP );
  poMacro( RIA );
  poMacro( LIA );
  poMacro( RSA );
  poMacro( LSA );
  poMacro( IRP );
  poMacro( ILP );
  poMacro( SRP );
  poMacro( SLP );
  poMacro( IRA );
  poMacro( ILA );
  poMacro( SRA );
  poMacro( SLA );
  poMacro( RPI );
  poMacro( LPI );
  poMacro( RAI );
  poMacro( LAI );
  poMacro( RPS );
  poMacro( LPS );
  poMacro( RAS );
  poMacro( LAS );
  poMacro( PRI );
  poMacro( PLI );
  poMacro( ARI );
  poMacro( ALI );
  poMacro( PRS );
  poMacro( PLS );
  poMacro( ARS );
  poMacro( ALS );
  poMacro( IPR );
  poMacro( SPR );
  poMacro( IAR );
  poMacro( SAR );
  poMacro( IPL );
  poMacro( SPL );
  poMacro( IAL );
  poMacro( SAL );
  poMacro( PIR );
  poMacro( PSR );
  poMacro( AIR );
  poMacro( ASR );
  poMacro( PIL );
  poMacro( PSL );
  poMacro( AIL );
  poMacro( ASL );
}

int main( int argc, char ** argv )
{
  std::string inputFileName;
  std::string outputFileName;
  std::string desiredString;
  std::string givenString;

  std::stringstream description;
  description << "Re-orient image." << std::endl;
  description << "Orientation is defined by a combination of three characters:" << std::endl;
  description << "\tL,R = left, right" << std::endl;
  description << "\tA,P = anterior, posterior" << std::endl;
  description << "\tS,I = superior, inferior" << std::endl;
  description << "e.g., AIL: x = A-P, y = I-S, z = L-R";

  tkd::CmdParser p( "orient", description.str() );

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetInput( "filename" )
    ->SetDescription( "Input image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetInput( "filename" )
    ->SetDescription( "Output image" )
    ->SetRequired( true );

  p.AddArgument( desiredString, "desired" )
    ->AddAlias( "d" )
    ->SetInput( "orientation" )
    ->SetDescription( "Desired output orientation (default: RAI)" );

  p.AddArgument( givenString, "given" )
    ->AddAlias( "g" )
    ->SetInput( "orientation" )
    ->SetDescription( "Given input orientation (default: determine from image)" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  itk::SpatialOrientation::ValidCoordinateOrientationFlags given, desired;
  given = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID;
  desired = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;

  if ( desiredString != "" )
    {
    ParseOrientation( desiredString, desired );
    }

  if ( givenString != "" )
    {
    ParseOrientation( givenString, given );
    }

  itk::ImageIOBase::Pointer io =
    itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(), itk::ImageIOFactory::ReadMode );

  if ( !io )
    {
    std::cerr << "Could not create a reader for " << inputFileName << std::endl;
    return -1;
    }

  io->SetFileName( inputFileName.c_str() );
  io->ReadImageInformation();

#define switchMacro( itkPixel, pixel, dimension ) \
  if ( io->GetNumberOfDimensions() == dimension && io->GetComponentType() == itk::ImageIOBase::itkPixel ) \
  { \
    return OrientImage< pixel, dimension >( inputFileName, outputFileName ); \
  }

//  switchMacro( FLOAT, float, 3 );
//  switchMacro( FLOAT, float, 4 );

  return OrientImage< float, 3 >( inputFileName, outputFileName, desired, given );

  std::cerr << "Unsupported pixel/dimension" << std::endl;
  return -1;
}
