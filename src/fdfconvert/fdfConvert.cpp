#include "fdfConvert.h"
#include "itkImageFileWriter.h"
#include "itksys/SystemTools.hxx"
#include "itksys/Directory.hxx"
#include "itksys/RegularExpression.hxx"
#include <fstream>
#include <iostream>

namespace fdf
{

FDFConvert::FDFConvert()
{
  m_MinZ = itk::NumericTraits< double >::max();
  m_MaxZ = itk::NumericTraits< double >::min();
  m_Verbose = false;
}

void FDFConvert::SetVerbose( bool verbose )
{
  m_Verbose = verbose;
}

bool FDFConvert::ReadLine( const std::string& line, Dictionary& dictionary )
{
  // header line
  itksys::RegularExpression re( "(int|float|char)\\ *(\\*)?([a-zA-Z0-9_]+)(\\[\\])?\\ *=\\ *(\"(.+)\"\\ *;?)?(\\{(.+)\\}\\ *;?)?([^;]*)" );

  // matrix/array of values separated by ,
  itksys::RegularExpression re2( "^([^,\\ ]+)[,\\ ]*(.*)$" );

  // string expression between quotes
  itksys::RegularExpression re3( "^\"(.+)\"$" );

  if ( re.find( line ) )
    {
    // int|float|char
    std::string type = re.match( 1 );

    // * denotes array
    bool isArray = re.match( 2 ) == "*";

    // key name
    std::string key = re.match( 3 );

    // matrix denoted by []
    bool isMatrix = re.match( 4 ) == "[]";

    // value (either in {}, in "", or plain value)
    std::string value = re.match( 6 ) != "" ? re.match( 6 ) : ( re.match( 8 ) != "" ? re.match( 8 ) : re.match( 9 ) );

    // unquote
    if ( value == "\"\"" )
      {
      value = "";
      }

    std::vector< std::string > values;

    if ( isMatrix )
      {
      // get matrix values separated by comma's and possibly quoted
      while ( re2.find( value ) )
        {
        if ( re3.find( re2.match( 1 ) ) )
          {
          values.push_back( re3.match( 1 ) );
          }
        else
          {
          values.push_back( re2.match( 1 ) );
          }

        value = re2.match( 2 );
        }
      }
    else
      {
      values.push_back( value );
      }

    dictionary[ key ] = values;
    return true;
    }

  return false;
}

void FDFConvert::Read( const std::string& path )
{
  typedef std::map< int, std::string > SliceMap;
  typedef std::map< int, SliceMap > ImageMap;
  typedef std::map< int, ImageMap > EchoMap;

  EchoMap echos;

  std::set< int > sliceCount, echoCount, imageCount;

  // get directory files
  itksys::Directory dir;
  dir.Load( path.c_str() );

  // valid files; we only support 2D FDF-files
  itksys::RegularExpression re( "slice([0-9]+)image([0-9]+)echo([0-9]+)\\.fdf" );

  for( unsigned long i = 0; i < dir.GetNumberOfFiles(); ++i )
    {
    std::string filename = dir.GetFile( i );

    if ( re.find( filename ) )
      {
      // valid filename, get full path
      std::string fullPath = path + std::string( "/" ) + filename;
      std::string file = itksys::SystemTools::CollapseFullPath( fullPath.c_str() );

      // get image number, echo number, slice number
      std::stringstream ss;
      ss << re.match( 1 ) << " " << re.match( 2 ) << " " << re.match( 3 );

      int slice, image, echo;
      ss >> slice >> image >> echo;

      --slice;
      --image;
      --echo;

      // store values in an orderer map
      echos[ echo ][ image ][ slice ] = file;

      sliceCount.insert( slice );
      imageCount.insert( image );
      echoCount.insert( echo );
      }
    else
      {
      if ( m_Verbose ) std::cout << "Skip " << filename << std::endl;
      }
    }

  this->m_NumberOfEchos = echoCount.size();
  this->m_NumberOfImages = imageCount.size();
  this->m_NumberOfSlices = sliceCount.size();

  // read files
  for( EchoMap::iterator i = echos.begin(); i != echos.end(); ++i )
    {
    int echo = i->first;
    for( ImageMap::iterator j = i->second.begin(); j != i->second.end(); ++j )
      {
      int image = j->first;
      for( SliceMap::iterator k = j->second.begin(); k != j->second.end(); ++k )
        {
        int slice = k->first;

        if ( m_Verbose ) std::cout << "Read " << slice << "/" << image << "/" << echo << std::endl;

        LoadSlice( k->second, slice, image + echo * this->m_NumberOfImages );
        }
      }
    }
}

void FDFConvert::Write( const std::string& output )
{
  if ( m_Verbose ) std::cout << "Write " << output << std::endl;

  typedef itk::ImageFileWriter< SeriesType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( output.c_str() );
  writer->SetInput( m_Data );
  writer->Update();
}

void FDFConvert::LoadSlice( const std::string& filename, int slice, int volume )
{
  std::cout << "Read file " << filename << std::endl;

  Dictionary dictionary;

  // open file
  std::ifstream in( filename.c_str() );
  bool first = true;

  while( !in.eof() )
    {
    std::string line;
    std::getline( in, line );

    // skip first line
    if ( first )
      {
      first = false;
      continue;
      }

    // as long as we receive valid header lines, continue
    if ( !ReadLine( line, dictionary ) )
      {
      break;
      }
    }

  // get byte size (should be 4 for float)
  int byteSize = dictionary.Get< int >( "bits" ) / 8;

  // slice matrix size
  int matrixSize = dictionary.Get< int >( "matrix", 0 ) * dictionary.Get< int >( "matrix", 1 );

  // update number of slices
  m_NumberOfSlices = dictionary.Get< int >( "slices", 0, m_NumberOfSlices );

  // image size from expected number of slices
  int imageSize = matrixSize * m_NumberOfSlices;

  if ( byteSize != 4 )
    {
    throw itk::ExceptionObject( "Only 4 bytes float images are supported" );
    }

  if ( !m_Data )
    {
    // initialize 4D data set

    m_Data = SeriesType::New();
    SeriesType::RegionType region;
    SeriesType::SizeType size;

    // matrix size
    size[ 0 ] = dictionary.Get< int >( "matrix", 0, 0 );
    size[ 1 ] = dictionary.Get< int >( "matrix", 1, 0 );
    size[ 2 ] = m_NumberOfSlices;
    size[ 3 ] = m_NumberOfImages * m_NumberOfEchos;
    region.SetSize( size );
    m_Data->SetRegions( region );

    // spacing and origin
    SeriesType::SpacingType spacing;
    SeriesType::PointType origin;
    SeriesType::DirectionType direction = m_Data->GetDirection();

    spacing[ 0 ] = ( dictionary.Get< double >( "span", 0, 0 ) * 10.0 ) / dictionary.Get< double >( "matrix", 0, 1 );
    spacing[ 1 ] = ( dictionary.Get< double >( "span", 1, 0 ) * 10.0 ) / dictionary.Get< double >( "matrix", 1, 1 );
    spacing[ 2 ] = dictionary.Get< double >( "roi", 2, 0 ) * 10.0;
    spacing[ 3 ] = 1;

    origin[ 0 ] = dictionary.Get< double >( "origin", 0, 0 ) + dictionary.Get< double >( "location", 0, 0 );
    origin[ 1 ] = dictionary.Get< double >( "origin", 1, 0 ) + dictionary.Get< double >( "location", 1, 0 );
    origin[ 2 ] = 0;
    origin[ 3 ] = 0;

    // direction
    for( int x = 0; x < 3; ++x )
      {
      for( int y = 0; y < 3; ++y )
        {
        direction[ x ][ y ] = dictionary.Get< double >( "orientation", x * 3 + y, 0 );
        }
      }

    m_Data->SetSpacing( spacing );
    m_Data->SetDirection( direction );
    m_Data->SetOrigin( origin );

    // allocate memory
    m_Data->Allocate();

    // fill with zeros
    m_Data->FillBuffer( 0 );
    }

  // unused, but keep track of each slice's location
  double location = dictionary.Get< double >( "location", 2, 0 ) * 10.0;

  m_MinZ = m_MinZ < location ? m_MinZ : location;
  m_MaxZ = m_MaxZ > location ? m_MaxZ : location;

  // set buffer pointer to position of current slice
  float* buffer = m_Data->GetPixelContainer()->GetBufferPointer() + imageSize * volume + matrixSize * slice;

  // move file reader pointer
  in.seekg( -1 * sizeof( float ) * matrixSize, std::ios::end );

  // read data
  in.read( reinterpret_cast< char* >( buffer ), sizeof( float ) * matrixSize );
}

} // end namespace fdf
