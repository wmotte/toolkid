/*
 * copyinfo.cpp
 *
 *  Created on: Jul 27, 2009
 *      Author: wim
 */
#include "itkImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "tkdCmdParser.h"

#include <iostream>
#include <fstream>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

/**
 * Convert SignallExpress Voltage.txt ascii output (EEG monitoring; National Instruments).
 * To nii.gz Image.
 */
class Txt2Nii
{

protected:
	/**
	 * Replace all sub-strings.
	 */
	inline std::string replaceAll( const std::string& s, const std::string& f, const std::string& r )
	{
		if ( s.empty() || f.empty() || f == r || s.find( f ) == std::string::npos )
		{
			return s;
		}
		std::ostringstream build_it;
		size_t i = 0;
		for ( size_t pos; ( pos = s.find( f, i ) ) != std::string::npos; )
		{
			build_it.write( &s[i], pos - i );
			build_it << r;
			i = pos + f.size();
		}
		if ( i != s.size() )
		{
			build_it.write( &s[i], s.size() - i );
		}
		return build_it.str();
	}

public:

	/**
	 * Run app.
	 */
	void run( const std::string& input, const std::string& output, float seconds )
	{
		process( input, output, seconds );
	}

	/**
	 * New.
	 */
	std::vector< std::vector< double > > getChannels2( const std::string& inputFileName )
	{
		typedef boost::tokenizer< boost::char_separator< char > > TokenizerType;

		std::vector< std::vector< double > > data;

		std::ifstream in;
		in.open( inputFileName.c_str() );

		if ( !in )
		{
			std::cerr << "Unable to open: " << inputFileName;
			exit( EXIT_FAILURE );
		}

		std::string line;

		// header ...
		boost::regex dataTag( "data:" );
		boost::smatch m;
		boost::sregex_token_iterator end;

		if( UseRowNames() )
		// bet line ...
		while ( getline( in, line ) )
		{
			if ( boost::regex_search( line, m, dataTag ) )
			{
				break;
			}
		}

		// continue with data ...
		while ( getline( in, line ) )
		{
			// replace comma's with dots...
			boost::regex commaDotTag( ",", boost::regex_constants::icase | boost::regex_constants::perl );
			std::string newLine  = boost::regex_replace( line, commaDotTag, "." );

			try
			{
				std::vector< double > line_data;

				TokenizerType tok( newLine, boost::char_separator< char >( "\t" ) );

				for ( TokenizerType::iterator id = tok.begin(); id != tok.end(); ++id )
				{
					double value = boost::lexical_cast< double >( *id );
					line_data.push_back( value );
				}
				data.push_back( line_data );
			}
			catch ( boost::bad_lexical_cast& e )
			{
				// If parsing fails, try ss method ...
				return getChannels( inputFileName );
			}
		}

		return data;
	}

	/**
	 * Get channels as vectors from given ascii-file.
	 *
	 * First read header and determine number of channels.
	 */
	std::vector< std::vector< double > > getChannels( const std::string& inputFileName )
	{

		std::ifstream in;
		in.open( inputFileName.c_str() );

		if ( !in )
		{
			std::cerr << "Unable to open: " << inputFileName;
			exit( EXIT_FAILURE );
		}

		// first read header, which should end with "data: "...
		bool inHeader = true;
		while ( !in.eof() && in.good() && inHeader )
		{
			std::string line;
			std::getline( in, line );

			int pos = line.find( "data:" );

			if ( pos == 0 )
			{
				inHeader = false;
				break;
			}
		}

		std::vector< std::vector< double > > data;

		while ( !in.eof() && in.good() )
		{
			std::string line;
			std::getline( in, line );

			// replace komma's by dots...
			std::string goodLine = replaceAll( line, ",", "." );

			std::stringstream ss;
			ss.precision( 6 );
			ss << goodLine;

			std::vector< double > line_data;
			bool valid = false;
			while ( !ss.eof() && ss.good() )
			{
				double x;
				ss >> x;
				if ( ss.good() )
				{
					line_data.push_back( x );
					valid = true;
				}
			}

			if ( valid )
			{
				data.push_back( line_data );
			}
		}
		in.close();

		return data;
	}

	/**
	 * Write 4D image with each voxel a channel over time...
	 */
	void process( const std::string& inputFileName, const std::string& outputFileName, const float seconds )
	{
		if ( seconds > 30 )
		{
			std::cerr << "Series are restricted to < 30 seconds (NIFTI-header restriction)!" << std::endl;
			exit( EXIT_FAILURE );
		}
		unsigned int milliseconds = seconds * 1000;

		const unsigned int Dimension = 4;
		typedef double PixelType;
		typedef itk::Image< PixelType, Dimension > ImageType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		ImageType::Pointer image = ImageType::New();

		//std::vector< std::vector< double > > channels = getChannels( inputFileName );
		std::vector< std::vector< double > > channels = getChannels2( inputFileName );

		if ( !channels.empty() )
		{

			// set image size and allocate data...
			ImageType::RegionType region;
			ImageType::SizeType size;

			unsigned int numberOfChannels = channels[0].size();
			unsigned int totalSize = channels.size();
			unsigned int timeSize = totalSize;

			size[0] = numberOfChannels;
			size[1] = ( channels.size() / milliseconds ) + 1;
			size[2] = 1;
			size[3] = milliseconds;

			region.SetSize( size );
			image -> SetRegions( region );
			image -> Allocate();
			image -> FillBuffer( 0 );

			ImageType::IndexType index;

			// for all channels...
			for ( unsigned int channelIndex = 0; channelIndex < numberOfChannels; channelIndex++ )
			{

				// Cut data in blocks of 30 sec...
				unsigned int blockCount = -1;

				// each channel...
				for ( unsigned int timeIndex = 0; timeIndex < timeSize; timeIndex++ )
				{

					// switch to new voxel...
					if ( ( timeIndex % milliseconds ) == 0 )
					{
						blockCount++;
					}

					index[0] = channelIndex;
					index[1] = blockCount;
					index[2] = 0;
					index[3] = ( timeIndex - ( blockCount * milliseconds ) );

					std::vector< double > points = channels[timeIndex];
					double pixelValue = points[channelIndex];

					image -> SetPixel( index, pixelValue );
				}
			}

			WriterType::Pointer writer = WriterType::New();
			writer -> SetFileName( outputFileName );
			writer -> SetInput( image );

			try
			{
				writer -> Update();
			} catch ( itk::ExceptionObject & exp )
			{
				std::cerr << "ERROR - could not write: " << inputFileName << std::endl;
				exit( EXIT_FAILURE );
			}
		} else
		{
			std::cerr << "ERROR - No channels in: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		}
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] )
{
	std::string input;
	std::string output;
	float seconds = 30.0;

	tkd::CmdParser parser( argv[0], "EEG - Convert 'SignallExpress' ascii-file to NIFTI 4D Image" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input file (ascii)" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output file (image)" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( seconds, "seconds" ) -> AddAlias( "s" ) -> SetInput( "<float>" ) -> SetDescription(
			"Seconds per time-series (default 30)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Txt2Nii txt2nii = Txt2Nii();
	txt2nii.run( input, output, seconds );
}

