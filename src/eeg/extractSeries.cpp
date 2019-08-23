/*
 * copyinfo.cpp
 *
 *  Created on: Jul 27, 2009
 *      Author: wim
 */
#include "itkImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include "tkdCmdParser.h"

#include <iostream>
#include <fstream>
#include <string>

/**
 * Extract power-points from one series (normally 30 sec; EEG).
 */
class ExtractSeries {

public:

	/**
	 * Run app.
	 */
	void run( const std::string& input, const unsigned int channel, const unsigned int rangeOffset, const unsigned int rangeSize ) {

		const unsigned int Dimension = 4;
		typedef double PixelType;
		typedef itk::Image< PixelType, Dimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;

		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( input );
		ImageType::Pointer image = reader -> GetOutput();
		reader -> Update();

		ImageType::SizeType size = image -> GetLargestPossibleRegion().GetSize();
		ImageType::IndexType index = image -> GetLargestPossibleRegion().GetIndex();

		//std::cout << "Index: " << index << std::endl; // [ 0 0  0  0 ]
		//std::cout << "Size: " << size << std::endl; // [ 2 73 1 30 ]

		unsigned int maxChannels = size[0];
		unsigned int maxRange = size[1];
		unsigned int maxTimeSeries = size[3];

		// checks...
		if ( channel >= maxChannels ) {
			std::cerr << "Max channel is: " << maxChannels - 1 << std::endl;
			exit( EXIT_FAILURE );
		} else if ( rangeOffset + rangeSize >= maxRange ) {
			std::cerr << "Max range is: " << maxRange - 1 << std::endl;
			exit( EXIT_FAILURE );
		}

		for ( unsigned int timePoint = 0; timePoint < maxTimeSeries; timePoint++ ) {

			for ( unsigned int range = rangeOffset; range < ( rangeOffset + rangeSize ); range++ ) {

				index[0] = channel;
				index[1] = range;
				index[2] = 0;
				index[3] = timePoint;

				std::cout << image -> GetPixel( index ) << ",";
			}
			std::cout << std::endl;
		}
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;
	int channel;
	int offset;
	int size;

	tkd::CmdParser parser( argv[0], "EEG - Print range of time-series as CSV" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription(
			"Input power spectrum file (Image)" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( channel, "channel" ) -> AddAlias( "c" ) -> SetInput( "<int>" ) -> SetDescription( "Channel to parse" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( offset, "offset" ) -> AddAlias( "o" ) -> SetInput( "<int>" ) -> SetDescription(
			"Range Offset (zero-based)" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( size, "size" ) -> AddAlias( "s" ) -> SetInput( "<int>" ) -> SetDescription(
				"Range Size" ) -> SetRequired( true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ExtractSeries extractSeries = ExtractSeries();
	extractSeries.run( input, channel, offset, size );
}

