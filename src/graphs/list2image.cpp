#include "tkdCmdParser.h"

#include "graphCommon.h"

/**
 * List -> 3D image.
 */
class List2Image
{

public:

	typedef double ValueType;

	/**
	 * Run.
	 */
	void run( const std::string& input, const std::string& output,
			  const std::string& reference, bool similarityData )
	{
		if ( !similarityData )
		{
			graph::Graph< ValueType >::VectorType vector;

			graph::Graph< ValueType >::GetVectorFromFile( input, vector );
			graph::Graph< ValueType >::Fill3DImageWithValues( reference, output, vector );
		}
		else
		{
			graph::Graph< ValueType >::VectorType data;
			graph::Graph< ValueType >::VectorType indices;
			
			GetSimilarityLabels( input, indices, data );
			
			graph::Graph< ValueType >::Fill3DImageWithValues( reference, output, indices, data );
		}
	}
	
protected:
	
	/**
	 * Get data and indices from similarity output file.
	 */
	void GetSimilarityLabels( const std::string& inputFileName, graph::Graph< ValueType >::VectorType& indices, 
			graph::Graph< ValueType >::VectorType& data )
	{
		std::ifstream in( inputFileName.c_str() );
		std::ifstream inLineCount( inputFileName.c_str() );
		unsigned int numberOfLines = std::count( std::istreambuf_iterator< char >( inLineCount ), std::istreambuf_iterator< char >(), '\n' );
		inLineCount.close();
		
		indices = graph::Graph< ValueType >::VectorType( numberOfLines );
		data = graph::Graph< ValueType >::VectorType( numberOfLines );
		
		if ( !in )
		{
			std::cerr << "*** ERROR ***: Unable to open file \"" << inputFileName << "\" for input.\n";
			exit( EXIT_FAILURE );
		}

		std::string line;	

		try
		{
			unsigned int i = 0;
			while ( getline( in, line ) )
			{
				// skip potential comment lines (given with '%').
				if( line.substr( 0, 1 ) == "%" )
					continue;
				if( line.empty() )
					continue;
				
				graph::Graph< ValueType >::TokType tok( line, boost::char_separator< char >( "\t" ) );

				graph::Graph< ValueType >::TokType::iterator id = tok.begin();
				ValueType index = boost::lexical_cast< ValueType >( *id );	
				id++;
				ValueType label = boost::lexical_cast< ValueType >( *id );
				
				indices.put( i, index );
				data.put( i, label );
				
				i++;
			}
		}
		catch ( std::exception& e )
		{
			std::cerr << "*** ERROR ***: could not parse input. "
				"Are you sure columns are separated with the '\\t' char?" << std::endl;

			exit( EXIT_FAILURE );
		}

		in.close();
	}
};

/**
 * Main hierarchy.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "list2image", "Insert list into 3D image." );

	std::string inputFileName;
	std::string outputFileName;
	std::string maskImageFileName;
	bool similarityData;


	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "filename" ) ->SetDescription( "Input text file (single column)" ) ->SetRequired(
			true )->SetMinMax( 1, 1 );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output 3D image file" ) ->SetRequired(
			true ) ->SetMinMax( 1, 1 );

	p.AddArgument( maskImageFileName, "mask" ) ->AddAlias( "m" ) ->SetInput( "filename" ) ->SetDescription( "Mask 3D image file (reference)" ) ->SetRequired(
			false ) ->SetMinMax( 1, 1 );
	
	p.AddArgument( similarityData, "similarity-data" ) ->AddAlias( "s" ) ->SetInput( "bool" ) ->SetDescription( "Input file contains similarity data" ) ->SetRequired(
			false ) ->SetMinMax( 1, 1 );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	List2Image list2Image;

	list2Image.run( inputFileName, outputFileName, maskImageFileName, similarityData );

	return EXIT_SUCCESS;
}



