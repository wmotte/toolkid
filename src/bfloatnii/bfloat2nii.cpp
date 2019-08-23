#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include "itkFlipImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileWriter.h"
#include "tkdCmdParser.h"

#include "itkCastImageFilter.h"

#define logMacro(x) \
  { \
  if ( m_Verbose ) \
    { \
      std::cout << x << std::endl; \
    } \
  }

class ExportBFloat
{
private:
	bool m_Verbose;

	typedef float PixelType;
	typedef short ShortPixelType;
	static const unsigned int ImageDimension = 3;
	static const unsigned int SeriesDimension = 4;
	typedef itk::Image< PixelType, ImageDimension > ImageType;
	typedef itk::Image< PixelType, SeriesDimension > SeriesType;
	typedef itk::Image< ShortPixelType, ImageDimension > ShortImageType;
	typedef itk::Image< ShortPixelType, SeriesDimension > ShortSeriesType;

	void ByteSwap( char* block, int size )
	{
		register int i = 0;
		register int j = size - 1;

		for ( ; i < j; ++i, --j )
		{
			std::swap( block[i], block[j] );
		}
	}

	/**
	 * Add bshort support (wim; 27-03-2010).
	 */
	ImageType::Pointer ImportImage( std::string inputFileName )
	{
		ImageType::SizeType size;
		ImageType::RegionType region;
		ShortImageType::RegionType shortRegion;

		bool bshort = false;

		if ( inputFileName.length() > 7 && inputFileName.substr( inputFileName.length() - 6, 6 ) == "bfloat" )
		{
			inputFileName = inputFileName.substr( 0, inputFileName.length() - 6 );
		} else if ( inputFileName.length() > 7 && inputFileName.substr( inputFileName.length() - 6, 6 ) == "bshort" )
		{
			inputFileName = inputFileName.substr( 0, inputFileName.length() - 6 );
			bshort = true;
		}

		std::ifstream hin( ( inputFileName + std::string( "hdr" ) ).c_str() );
		if ( !hin.is_open() )
		{
			std::cerr << "Error reading header file '" << inputFileName << "hdr'" << std::endl;
			return 0;
		}

		size.Fill( 0 );
		for ( int i = 0; i < 3; ++i ) // TODO bshort/hdr Mandeville contains other header format! (4D)
		{
			if ( !hin.eof() )
			{
				hin >> size[i];
			}
		}

		hin.close();

		int tmp = size[1];
		size[1] = size[0];
		size[0] = tmp;

		std::ifstream in( ( inputFileName + std::string( bshort ? "bshort" : "bfloat" ) ).c_str(), std::ios_base::binary
				| std::ios_base::ate );

		if ( !in.is_open() )
		{
			std::cerr << "Error reading file '" << inputFileName << "'" << std::endl;
			return 0;
		}

		std::ifstream::pos_type streamSize = in.tellg();
		in.seekg( 0, std::ios_base::beg );

		ImageType::Pointer image = ImageType::New();

		if ( !bshort )
		{
			int bufferSize = streamSize / 4;

			float* data = new float[bufferSize];
			in.read( reinterpret_cast< char* > ( data ), streamSize );
			in.close();

			logMacro( "Byte swap block of size " << size << " bytes" );

			for ( int i = 0; i < bufferSize; ++i )
			{
				ByteSwap( reinterpret_cast< char* > ( &( data[i] ) ), 4 );
			}

			region.SetSize( size );
			image->SetRegions( region );

			if ( !image->GetPixelContainer() )
			{
				std::cout << "Warning: container is null" << std::endl;
			}

			image->GetPixelContainer()->SetImportPointer( data, bufferSize, true );

		} else // bshort
		{
			ShortImageType::Pointer shortImage = ShortImageType::New();

			int bufferSize = streamSize / 2;

			short* data = new short[bufferSize];
			in.read( reinterpret_cast< char* > ( data ), streamSize );
			in.close();

			logMacro( "Byte swap block of size " << size << " bytes" );

			for ( int i = 0; i < bufferSize; ++i )
			{
				//ByteSwap( reinterpret_cast< char* > ( &( data[i] ) ), 2 ); // -> Type 1... Mandeville
			}

			shortRegion.SetSize( size );
			shortImage->SetRegions( shortRegion );
			shortImage->GetPixelContainer()->SetImportPointer( data, bufferSize, true );

			typedef itk::CastImageFilter< ShortImageType, ImageType > CastImageFilterType;
			CastImageFilterType::Pointer convertFilter = CastImageFilterType::New();
			convertFilter->SetInput( shortImage );
			convertFilter->Update();
			image = convertFilter->GetOutput();
		}

		return image;
	}

	SeriesType::Pointer JoinImages( const std::vector< ImageType::Pointer >& images, bool seriesBySlices = false )
	{
		int numberOfImages = images.size();
		SeriesType::Pointer series = SeriesType::New();

		ImageType::RegionType imageRegion = images[0]->GetLargestPossibleRegion();
		ImageType::SizeType imageSize = imageRegion.GetSize();

		SeriesType::SizeType size;

		if ( seriesBySlices )
		{
			size[0] = imageSize[0];
			size[1] = imageSize[1];
			size[2] = numberOfImages;
			size[3] = imageSize[2];
		}
		else
		{
			size[0] = imageSize[0];
			size[1] = imageSize[1];
			size[2] = imageSize[2];
			size[3] = numberOfImages;
		}

		SeriesType::RegionType region;
		region.SetSize( size );

		series->SetRegions( region );
		series->Allocate();

		for ( int i = 0; i < numberOfImages; ++i )
		{
			SeriesType::RegionType seriesRegion = region;
			SeriesType::IndexType seriesIndex = seriesRegion.GetIndex();
			SeriesType::SizeType seriesSize = seriesRegion.GetSize();

			if ( seriesBySlices )
			{
				seriesIndex[2] = i;
				seriesSize[2] = 1;
			} else
			{
				seriesIndex[3] = i;
				seriesSize[3] = 1;
			}

			seriesRegion.SetSize( seriesSize );
			seriesRegion.SetIndex( seriesIndex );

			itk::ImageRegionConstIterator< ImageType > itImage( images[i], images[i]->GetLargestPossibleRegion() );
			itk::ImageRegionIterator< SeriesType > itSeries( series, seriesRegion );

			for ( ; !itImage.IsAtEnd(); ++itImage, ++itSeries )
			{
				itSeries.Set( itImage.Get() );
			}
		}

		return series;
	}

	bool ParseProcPar( std::string filename, std::map< std::string, std::vector< std::string > >& procpar )
	{
		if ( filename.length() < 7 || filename.substr( filename.length() - 7, 7 ) != "procpar" )
		{
			if ( filename.length() < 4 || ( filename.substr( filename.length() - 4, 4 ) != ".fid" && filename.substr(
					filename.length() - 4, 4 ) != "fid/" ) )
			{
				std::stringstream ss;
				ss << filename << ".fid/procpar";
				filename = ss.str();
			} else
			{
				std::string separator = "/";
				if ( filename.substr( filename.length() - 1, 1 ) == "/" )
				{
					separator = "";
				}

				std::stringstream ss;
				ss << filename << separator << "procpar";
				filename = ss.str();
			}
		}

		logMacro( "Read procpar file from " << filename );

		std::ifstream in( filename.c_str() );
		if ( !in.is_open() )
		{
			std::cerr << "Could not read " << filename << std::endl;
			return false;
		}

		while ( !in.eof() )
		{
			std::string id;
			std::string stuff;

			in >> id;

			if ( id == "" )
			{
				break;
			}

			for ( int i = 0; i < 10; ++i )
			{
				in >> stuff;
			}

			int numberOfValues = 0;

			for ( int j = 0; j < 2; ++j )
			{
				in >> numberOfValues;

				for ( int i = 0; i < numberOfValues; ++i )
				{
					std::string value;
					in >> value;

					if ( value[0] == '"' )
					{
						int delimiters;

						do
						{
							delimiters = 0;

							for ( unsigned int k = 0; k < value.length(); ++k )
							{
								if ( value[k] == '"' )
								{
									++delimiters;
								}
							}

							if ( delimiters % 2 != 0 )
							{
								std::string add;
								in >> add;
								value = value + std::string( " " ) + add;
							}
						} while ( delimiters % 2 != 0 && !in.eof() );
					}

					logMacro( "Procpar '" << id << "' (" << ( i + 1 ) << "/" << numberOfValues << ") = " << value );

					procpar[id].push_back( value );
				}
			}
		}

		in.close();

		return true;
	}

	ImageType::Pointer FlipImage( ImageType::Pointer input, bool flipX, bool flipY, bool flipZ, bool flipOrigin = false )
	{
		typedef itk::FlipImageFilter< ImageType > FlipType;
		FlipType::Pointer flip = FlipType::New();
		flip->SetInput( input );
		FlipType::FlipAxesArrayType f;
		f[0] = flipX;
		f[1] = flipY;
		f[2] = flipZ;
		flip->SetFlipAxes( f );
		flip->SetFlipAboutOrigin( flipOrigin );

		typedef itk::ChangeInformationImageFilter< ImageType > ChangeType;
		ChangeType::Pointer change = ChangeType::New();
		change->SetInput( flip->GetOutput() );
		change->ChangeOriginOn();
		change->ChangeSpacingOn();
		change->ChangeDirectionOn();

		ImageType::PointType origin;
		ImageType::SpacingType spacing;
		ImageType::DirectionType direction;

		origin.Fill( 0 );
		spacing.Fill( 1 );

		for ( int i = 0; i < 3; ++i )
		{
			for ( int j = 0; j < 3; ++j )
			{
				direction[i][j] = ( i == j ? 1 : 0 );
			}
		}

		change->SetOutputOrigin( origin );
		change->SetOutputSpacing( spacing );
		change->SetOutputDirection( direction );
		change->Update();

		return change->GetOutput();
	}

public:

	int Run( int argc, char ** argv )
	{
		tkd::CmdParser p( "bfloat2nii", "Import bfloat images to nifti" );

		std::vector< std::string > filenames;
		std::string outputFileName;
		bool seriesBySlices = false;
		m_Verbose = false;
		bool flipY = false;
		std::string procparPath;
		std::vector< double > inputSpacing;
		std::vector< double > inputOrigin;
		std::vector< double > scaleSpacing;

		p.AddArgument( filenames, "input" ) ->AddAlias( "i" ) ->SetInput( "filename(s)" ) ->SetDescription( "Input bfloat images" ) ->SetRequired(
				true );

		p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output nifti image" ) ->SetRequired(
				true );

		p.AddArgument( seriesBySlices, "slice-series" ) ->AddAlias( "s" ) ->SetInput( "boolean" ) ->SetDescription(
				"Input images are 3D time series of a single slice" );

		p.AddArgument( m_Verbose, "verbose" ) ->AddAlias( "v" ) ->SetInput( "boolean" ) ->SetDescription( "Verbose output" );

		p.AddArgument( flipY, "flip-y" ) ->AddAlias( "fy" ) ->SetInput( "boolean" ) ->SetDescription( "Flip y-axis" );

		p.AddArgument( procparPath, "procpar" ) ->AddAlias( "p" ) ->SetInput( "path" ) ->SetDescription( "Path to procpar" );

		p.AddArgument( inputSpacing, "spacing" ) ->AddAlias( "sp" ) ->SetInput( "x y z [t]" ) ->SetDescription( "Image voxel spacing" ) ->SetMinMax(
				3, 4 );

		p.AddArgument( inputOrigin, "origin" ) ->AddAlias( "or" ) ->SetInput( "x y z [t]" ) ->SetDescription( "Image origin" ) ->SetMinMax(
				3, 4 );

		p.AddArgument( scaleSpacing, "scale-spacing" ) ->AddAlias( "ss" ) ->SetInput( "x y z [t]" ) ->SetDescription(
				"Scale image voxel spacing" ) ->SetMinMax( 3, 4 );

		if ( !p.Parse( argc, argv ) )
		{
			p.PrintUsage( std::cout );
			return -1;
		}

		std::vector< ImageType::Pointer > images;
		for ( unsigned int i = 0; i < filenames.size(); ++i )
		{
			logMacro( "Read " << filenames[ i ] );
			ImageType::Pointer image = ImportImage( filenames[i] );
			if ( flipY )
			{
				image = FlipImage( image, false, true, false );
				if ( !image )
				{
					return -1;
				}
			}
			images.push_back( image );
		}

		logMacro( "Join 4D" );
		SeriesType::Pointer series = JoinImages( images, seriesBySlices );
		SeriesType::SpacingType spacing = series->GetSpacing();
		SeriesType::PointType origin = series->GetOrigin();
		SeriesType::SizeType size = series->GetLargestPossibleRegion().GetSize();

		if ( procparPath != "" )
		{
			std::map< std::string, std::vector< std::string > > procpar;
			ParseProcPar( procparPath, procpar );

			std::string sequence = procpar["seqfil"].front();
			sequence = sequence.substr( 1, sequence.length() - 2 );

			logMacro( "Sequence: " << sequence );

			float lro, lpe, slcdist;

			if ( sequence == "fse3d_ir" || sequence == "fse3d" || sequence == "gemsangio3Dnew" || sequence == "ge3d" )
			{
				lro = atof( procpar["lro"].front().c_str() );
				lpe = atof( procpar["lpe"].front().c_str() );
				float lpe2 = atof( procpar["lpe2"].front().c_str() );
				slcdist = lpe2 / static_cast< float > ( size[seriesBySlices ? 3 : 2] );
				logMacro( "lpe2 = " << lpe2 );
			} else
			{
				lro = atof( procpar["lro"].front().c_str() );
				lpe = atof( procpar["lpe"].front().c_str() );
				std::vector< std::string > pss = procpar["pss"];
				std::vector< float > spacings;
				for ( unsigned int i = 0; i < pss.size(); ++i )
				{
					spacings.push_back( atof( pss[i].c_str() ) );
					logMacro( "pss[" << i << "] = " << pss[ i ] );
				}

				sort( spacings.begin(), spacings.end() );

				float sum = 0;

				for ( unsigned int i = 1; i < pss.size(); ++i )
				{
					sum += ( spacings[i] - spacings[i - 1] );
				}

				slcdist = sum / static_cast< float > ( pss.size() - 1 );
			}

			spacing[0] = lro / static_cast< float > ( size[0] ) * 10.0;
			spacing[1] = lpe / static_cast< float > ( size[1] ) * 10.0;
			spacing[2] = slcdist * 10.0;

			logMacro( "lro = " << lro );
			logMacro( "lpe = " << lpe );
			logMacro( "slcdist = " << slcdist );

		}

		spacing[3] = 1;

		for ( unsigned int i = 0; i < 4; ++i )
		{
			spacing[i] = ( i < inputSpacing.size() ? inputSpacing[i] : spacing[i] ) * ( i < scaleSpacing.size() ? scaleSpacing[i] : 1.0 );
			origin[i] = ( i < inputOrigin.size() ? inputOrigin[i] : origin[i] );
		}

		logMacro( "Spacing: " << spacing );
		logMacro( "Origin: " << origin );

		typedef itk::ChangeInformationImageFilter< SeriesType > ChangeType;
		ChangeType::Pointer change = ChangeType::New();
		change->SetInput( series );
		change->ChangeSpacingOn();
		change->ChangeOriginOn();
		change->SetOutputSpacing( spacing );
		change->SetOutputOrigin( origin );

		typedef itk::ImageFileWriter< SeriesType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput( change->GetOutput() );
		writer->SetFileName( outputFileName.c_str() );
		writer->Update();

		return 0;
	}

};

int main( int argc, char ** argv )
{
	ExportBFloat e;

	return e.Run( argc, argv );
}

