#include "fidCommon.h"
#include "itksys/SystemTools.hxx"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include <fstream>
#include <map>

namespace fid
{

std::string Common::GetAbsolutePath( const std::string& execPath, const std::string& fileName )
{
#ifdef _WIN32
  std::string directorySeparator = "\\";
#else
  std::string directorySeparator = "/";
#endif

  std::string pathOut;
  std::string errorMsg;

  if ( itksys::SystemTools::FindProgramPath( execPath.c_str(), pathOut, errorMsg ) )
    {
    std::string programPath = itksys::SystemTools::GetProgramPath( pathOut.c_str() );
    programPath = itksys::SystemTools::GetParentDirectory( programPath.c_str() );
    std::string fullPath = programPath + directorySeparator + fileName;
    return fullPath;
    }

  return fileName;
}

bool Common::LoadTable( const std::string& fileName, vnl_matrix< int >& table )
{
  std::ifstream tableIn( fileName.c_str() );

  if ( !tableIn.good() )
    {
    return false;
    }

  int np, nv;
  tableIn >> np >> nv;

  if ( tableIn.fail() || tableIn.bad() || tableIn.eof() )
    {
    return false;
    }

  table.set_size( np, nv );

  for( int i = 0; i < np; ++i )
    {
    for( int j = 0; j < nv; ++j )
      {
      double x;
      tableIn >> x;
      if ( !tableIn.good() || tableIn.fail() || tableIn.bad() )
        {
        return false;
        }
      table( i, j ) = static_cast< int >( x ) - 1;
      }
    }

  tableIn.close();

  return true;
}

void Common::SaveImage( ImageType::Pointer image, const std::string& filename )
{
  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( filename.c_str() );
  writer->Update();
}

void Common::SaveSeries( SeriesType::Pointer series, const std::string& filename )
{
  SeriesType::SizeType size = series->GetLargestPossibleRegion().GetSize();
  if ( size[ 3 ] < 2 )
    {
    SaveImage( series, filename, 0 );
    return;
    }

  typedef itk::ImageFileWriter< SeriesType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( series );
  writer->SetFileName( filename.c_str() );
  writer->Update();
}

void Common::SaveImage( SeriesType::Pointer image, const std::string& filename, int imageIndex )
{
  typedef itk::ExtractImageFilter< SeriesType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  SeriesType::RegionType region = image->GetLargestPossibleRegion();
  SeriesType::IndexType index = region.GetIndex();
  SeriesType::SizeType size = region.GetSize();

  size[ 3 ] = 0;
  index[ 3 ] = imageIndex;

  region.SetSize( size );
  region.SetIndex( index );

  filter->SetInput( image );
  filter->SetExtractionRegion( region );
  filter->Update();

  SaveImage( filter->GetOutput(), filename );
}

Common::SeriesType::Pointer Common::ComplexToMagnitude( ComplexSeriesType::Pointer input )
{
  SeriesType::Pointer output = SeriesType::New();
  SeriesType::RegionType region = input->GetLargestPossibleRegion();
  output->CopyInformation( input );
  output->SetRegions( region );
  output->Allocate();

  itk::ImageRegionIterator< SeriesType > itOut( output, region );
  itk::ImageRegionConstIterator< ComplexSeriesType > itIn( input, region );

  for( ; !itOut.IsAtEnd(); ++itOut, ++itIn )
    {
    const ComplexType& pixel = itIn.Value();
    itOut.Set( vcl_sqrt( pixel.real() * pixel.real() + pixel.imag() * pixel.imag() ) );
    }

  return output;
}

Common::SeriesType::Pointer Common::ComplexToReal( ComplexSeriesType::Pointer input )
{
  SeriesType::Pointer output = SeriesType::New();
  SeriesType::RegionType region = input->GetLargestPossibleRegion();
  output->CopyInformation( input );
  output->SetRegions( region );
  output->Allocate();

  itk::ImageRegionIterator< SeriesType > itOut( output, region );
  itk::ImageRegionConstIterator< ComplexSeriesType > itIn( input, region );

  for( ; !itOut.IsAtEnd(); ++itOut, ++itIn )
    {
    itOut.Set( itIn.Value().real() );
    }

  return output;

}

Common::SeriesType::Pointer Common::ComplexToImaginary( ComplexSeriesType::Pointer input )
{
  SeriesType::Pointer output = SeriesType::New();
  SeriesType::RegionType region = input->GetLargestPossibleRegion();
  output->CopyInformation( input );
  output->SetRegions( region );
  output->Allocate();

  itk::ImageRegionIterator< SeriesType > itOut( output, region );
  itk::ImageRegionConstIterator< ComplexSeriesType > itIn( input, region );

  for( ; !itOut.IsAtEnd(); ++itOut, ++itIn )
    {
    itOut.Set( itIn.Value().imag() );
    }

  return output;
}

std::vector< float > Common::GetSlices( const Procparser& pp )
{
  int numberOfZ = pp.GetSize( "pss" );
  std::vector< float > slices;
  for( int i = 0; i < numberOfZ; ++i )
    {
    slices.push_back( pp.GetAs< float >( "pss", i ) );
    }

  return slices;
}

std::vector< int > Common::GetSliceTable( const std::vector< float >& pss )
{
  std::map< float, int > sliceOrder;
  for( unsigned int i = 0; i < pss.size(); ++i )
    {
    sliceOrder[ pss[ i ] ] = i;
    }

  std::vector< int > slices;
  for( std::map< float, int >::iterator i = sliceOrder.begin(); i != sliceOrder.end(); ++i )
    {
    slices.push_back( i->second );
    }

  return slices;
}

Common::ComplexSeriesType::Pointer Common::ReorientCoronal( ComplexSeriesType::Pointer input, const Procparser& pp, bool correctGradientDirection )
{
  /**
   * Output = RIP:
   * X+ -> RL
   * Y+ -> IS
   * Z+ -> PA
   */
  std::string orient = pp.GetAs< std::string >( "orient" );

  // find out whether read-out or phase-encoding gradients run in opposite direction
  if ( correctGradientDirection )
    {
    double strengthReadOutGradient = pp.GetAs< double >( "gro" );
    double strengthPhaseEncodingGradient = pp.GetAs< double >( "gpe" );

    if ( strengthReadOutGradient < 0 || strengthPhaseEncodingGradient < 0 )
      {
      bool flip[] = { strengthPhaseEncodingGradient < 0, strengthReadOutGradient < 0, false, false };
      input = fid::Common::Flip< ComplexType, 4 >( input, flip, false );
      }
    }

  if ( orient == "cor" )
    {
    bool flip[] = { true, false, false, false };
    int permutation[] = { 0, 2, 1, 3 };

    input = fid::Common::Permute< ComplexType, 4 >(
        fid::Common::Flip< ComplexType, 4 >( input, flip, false ),
        permutation );
    }
  else if ( orient == "cor90" )
    {
    bool flip[] = { false, false, false, false };
    int permutation[] = { 1, 2, 0, 3 };

    input = fid::Common::Permute< ComplexType, 4 >(
        fid::Common::Flip< ComplexType, 4 >( input, flip, false ),
        permutation );
    }
  else if ( orient == "trans90" )
    {
    bool flip[] = { true, false, false, false };
    int permutation[] = { 1, 0, 2, 3 };
    input = fid::Common::Permute< ComplexType, 4 >(
        fid::Common::Flip< ComplexType, 4 >( input, flip, false ),
        permutation );
    }
  else if ( orient == "trans" )
    {
    bool flip[] = { false, true, false, false };
    input = fid::Common::Flip< ComplexType, 4 >( input, flip, false );
    }
  else if ( orient == "sag" )
    {
    bool flip[] = { true, false, false, false };
    int permutation[] = { 2, 0, 1, 3 };
    input = fid::Common::Permute< ComplexType, 4 >(
        fid::Common::Flip< ComplexType, 4 >( input, flip, false ),
        permutation );
    }
  else
    {
    std::cout << "WARNING: unsupported orientation: " << orient << std::endl;
    }

  ComplexSeriesType::DirectionType direction;
  direction.SetIdentity();

//  direction( 0, 0 ) = 1;
//  direction( 0, 1 ) = 0;
//  direction( 0, 2 ) = 0;
//  direction( 0, 3 ) = 0;
//  direction( 1, 0 ) = 0;
//  direction( 1, 1 ) = 0;
//  direction( 1, 2 ) = -1;
//  direction( 1, 3 ) = 0;
//  direction( 2, 0 ) = 0;
//  direction( 2, 1 ) = 1;
//  direction( 2, 2 ) = 0;
//  direction( 2, 3 ) = 0;
//  direction( 3, 0 ) = 0;
//  direction( 3, 1 ) = 0;
//  direction( 3, 2 ) = 0;
//  direction( 3, 3 ) = 1;

  input->SetDirection( direction );

  return input;
}

} // end namespace fid
