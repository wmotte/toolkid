#include "graphCommon.h"

#include <sstream>
#include <map>
#include <utility> // make_pair
#include <boost/program_options.hpp>

#include <boost/tokenizer.hpp>

#include <boost/lexical_cast.hpp>

/**
 * VisabilityGrapMap.
 */
class VGMap
{

public:

	typedef double ValueType;

	/**
	 * Convert all mean label signal to visability matrices.
	 */
	void run( const std::string& timeSeriesFileName, const std::string& outputRoot, const std::string& labelFileName,
			std::vector< int > labelRange, bool weighted, std::vector< unsigned int > ignore, bool verbose, const std::string& labelNames )
	{

		// [ 1 ]: dimension checks ...
		unsigned int dims = graph::Graph< ValueType >::GetImageDimensions( timeSeriesFileName );

		if ( dims != 4 )
		{
			std::cerr << "*** ERROR ***: File " << timeSeriesFileName << " is not 4D!\n" << std::endl;
			exit( EXIT_FAILURE );
		}
		dims = graph::Graph< ValueType >::GetImageDimensions( labelFileName );

		if ( dims != 3 )
		{
			std::cerr << "*** ERROR ***: File " << labelFileName << " is not 3D!\n" << std::endl;
			exit( EXIT_FAILURE );
		}

		// [ 2 ]: create visability matrix for each label ...
		process( timeSeriesFileName, outputRoot, labelFileName, labelRange, weighted, ignore, verbose, labelNames );
	}

protected:

	/**
	 *
	 */
	void process( const std::string& timeSeriesFileName, const std::string& outputRoot, const std::string& labelFileName,
			std::vector< int > labelRange, bool weighted, std::vector< unsigned int > ignore, bool verbose, const std::string& labelNames )
	{

		typedef itk::Image< ValueType, 4 > InputType;
		typedef itk::Image< ValueType, 4 > OutputType;
		typedef itk::Image< int, 3 > LabelType;

		typedef itk::ImageFileReader< InputType > InputReaderType;
		typedef itk::ImageFileReader< LabelType > LabelReaderType;
		typedef itk::ImageFileWriter< OutputType > WriterType;

		InputReaderType::Pointer inputReader = InputReaderType::New();
		inputReader->SetFileName( timeSeriesFileName.c_str() );
		inputReader->Update();

		InputType::Pointer input = inputReader->GetOutput();
		inputReader = 0;

		LabelReaderType::Pointer labelReader = LabelReaderType::New();
		labelReader->SetFileName( labelFileName.c_str() );
		labelReader->Update();

		LabelType::Pointer label = labelReader->GetOutput();
		labelReader = 0;

		std::map< int, int > labels;

		if ( labelRange.size() == 2 && labelRange[1] >= labelRange[0] )
		{
			int count = 0;
			for ( int i = labelRange[0]; i <= labelRange[1]; i++ )
			{
				labels[i] = count++;
			}
		} else
		{
			int count = 0;
			itk::ImageRegionIterator< LabelType > it( label, label->GetLargestPossibleRegion() );
			while ( !it.IsAtEnd() )
			{
				if ( it.Value() > 0 && labels.find( it.Value() ) == labels.end() )
				{
					labels[it.Value()] = count++;
				}
				++it;
			}
		}

		int numberOfLabels = labels.size();

		InputType::RegionType region = input->GetLargestPossibleRegion();
		InputType::SizeType size = region.GetSize();
		int voxels = size[0] * size[1] * size[2];
		int timepoints = size[3];

		vnl_matrix_ref< ValueType > matrix( timepoints, voxels, input->GetPixelContainer()->GetBufferPointer() );
		const int* mask = label->GetPixelContainer()->GetBufferPointer();

		vnl_matrix< ValueType > mean( timepoints, numberOfLabels );
		mean.fill( 0 );

		vnl_vector< int > count( numberOfLabels );

		count.fill( 0 );

		for ( int i = 0; i < voxels; ++i )
		{
			if ( mask[i] < 1 )
			{
				continue;
			}

			int column = labels[mask[i]];

			mean.set_column( column, mean.get_column( column ) + matrix.get_column( i ) );
			++count[column];
		}

		for ( int i = 0; i < numberOfLabels; ++i )
		{
			mean.set_column( i, mean.get_column( i ) / static_cast< ValueType > ( count( i ) ) );
		}

		// extract ignored samples ...
		int vectorLength = timepoints;
		int vectorStart = 0;

		if ( ignore.size() == 2 )
		{
			vectorStart = ignore[0];
			vectorLength = timepoints - vectorStart - ignore[1];
		}

		mean = mean.extract( vectorLength, numberOfLabels, vectorStart, 0 );

		graph::Graph< ValueType >::MapType labelMap;
		graph::Graph< ValueType >::GetLabelNames( labelMap, labelNames, "\t" );

		if ( verbose )
			std::cout << "Constructing visability graph for label: " << std::endl;

		for ( int i = 1; i <= numberOfLabels; i++ )
		{
			vnl_vector< ValueType > labelTimeSeries = mean.get_column( i );

			std::string name;
			graph::Graph< ValueType >::GetLabelName( outputRoot, ".nii.gz", labelMap, name, i );

			if ( verbose )
			{
				std::cout << i;
				std::cout.flush();
				std::cout << '\r';
			}

			graph::VisabilityGraph< ValueType >::WriteVisabilityMatrix( labelTimeSeries, name, false, weighted );
		}
	}
};

// ************************************************************************************

/**
 * Throw error if required option is not specified.
 */
void required_option( const boost::program_options::variables_map& vm, const std::string& required_option )
{
	if ( vm.count( required_option ) == 0 )
		throw std::logic_error( "Option: '" + required_option + "' is required!" );
}

/**
 * Visability graph mapping using labels.
 */
int main( int argc, char* argv[] )
{
	namespace po = boost::program_options;

	// application description ...
	std::string description = "Create visability graphs from 4D time-series using labels.\n";

	// options ...
	std::string input;
	std::string output;
	std::string labels;
	std::string labelNames;

	std::vector< int > labelRange;
	std::vector< unsigned int > ignore;

	bool weighted;
	bool verbose;

	try
	{
		po::options_description desc( "Available options" );

		desc.add_options()

		( "input,i", po::value< std::string >( &input ), "string: 4D image time-series input file." )

		( "output,o", po::value< std::string >( &output ), "string: output root." )

		( "labels,l", po::value< std::string >( &labels ), "string: 3D labels file." )

		( "label-names,n", po::value< std::string >( &labelNames ),
				"string: 3D labels file with id and name (column format: e.g. 3 S1FL le)." )

		( "weighted,w", po::value< bool >( &weighted ) ->default_value( false ) ->zero_tokens(),
				"bool: use weighted visability graph construction (slope)." )

		( "verbose,v", po::value< bool >( &verbose ) ->default_value( false ) ->zero_tokens(), "bool: print progress." )

		( "ignore,d", po::value< std::vector< unsigned int > >( &ignore )->multitoken(), "vector: samples to ignore at begin and end." )

		( "label-range,r", po::value< std::vector< int > >( &labelRange )->multitoken(), "vector: specify a label range." )

		( "help,h", "bool: produce help message." );

		po::variables_map vm;
		po::store( po::parse_command_line( argc, argv, desc ), vm );
		po::notify( vm );

		// help message ...
		if ( vm.count( "help" ) )
		{
			std::cout << argv[0] << ": " << description << std::endl;
			std::cout << desc << "\n";

			return EXIT_SUCCESS;
		}

		// required options ...
		required_option( vm, "input" );
		required_option( vm, "output" );
		required_option( vm, "labels" );

		// run application ...
		VGMap vgMap;

		vgMap.run( input, output, labels, labelRange, weighted, ignore, verbose, labelNames );

	} catch ( std::exception& e )
	{
		std::cerr << "*** ERROR ***: " << e.what() << "\n";
		std::cerr << "Use \"" << argv[0] << " --help\" for information about application usage." << std::endl;

		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
