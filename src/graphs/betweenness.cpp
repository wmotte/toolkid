#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "itkNumericTraits.h"

#include "tkdCmdParser.h"

#include <iostream>
#include <utility>
#include <fstream>

#include <boost/utility.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

/**
 * (Relative) Betweenness centrality.
 */
template< class DirectedProperty = boost::undirectedS >
class BC
{
public:

	typedef double PixelType;

	typedef boost::property< boost::vertex_index_t, unsigned int, boost::property< boost::vertex_centrality_t, PixelType > > VertexPropertyType;
	typedef boost::property< boost::edge_weight_t, PixelType, boost::property< boost::edge_centrality_t, PixelType > > EdgePropertyType;
	typedef boost::adjacency_matrix< DirectedProperty, VertexPropertyType, EdgePropertyType > GraphMatrixType;

	typedef itk::Image< PixelType, 2 > ImageType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef vnl_vector< PixelType > VectorType;

	/**
	 * Run.
	 */
	void Run( const std::string& filename, PixelType lowerThreshold, PixelType upperThreshold, const std::string& outputFileName,
			bool relative, bool dominance )
	{
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		typename ImageType::Pointer image = reader->GetOutput();
		reader = 0;

		PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
		typename ImageType::RegionType region = image->GetLargestPossibleRegion();
		typename ImageType::SizeType size = region.GetSize();
		int rows = size[0];
		int cols = size[1];

		DataMatrixType data( rows, cols, buffer );

		typename ImageType::Pointer output = ImageType::New();
		output->CopyInformation( image );
		output->SetRegions( image->GetLargestPossibleRegion() );
		output->Allocate();
		output->FillBuffer( 0 );

		DataMatrixType distance( rows, cols, output->GetPixelContainer()->GetBufferPointer() );

		GraphMatrixType g( rows );

		/* TEST List -> weighted version works, and indeed lower values are shorter paths,
		GraphMatrixType g( 5 );
		boost::add_edge( 0, 1, EdgePropertyType( 1.0 ), g );
		boost::add_edge( 1, 4, EdgePropertyType( 0.2 ), g );
		boost::add_edge( 4, 3, EdgePropertyType( 0.3 ), g );
		boost::add_edge( 0, 2, EdgePropertyType( 1.0 ), g );
		boost::add_edge( 2, 3, EdgePropertyType( 1.0 ), g );
		*/

		const bool isUndirected = typeid(DirectedProperty) == typeid(boost::undirectedS);
		for ( int i = 0; i < rows; ++i )
		{
			int start = isUndirected ? i + 1 : 0;
			for ( int j = start; j < cols; ++j )
			{
				if ( i == j )
				{
					continue;
				}

				PixelType weight = data( i, j );

				if ( weight <= upperThreshold && weight > 0 && weight > lowerThreshold )
				{
					// bc makes use of fwsp, so invert weight (correlation coefficients ...)
					boost::add_edge( i, j, EdgePropertyType( 1.0 / weight ), g );
				}
			}
		}

		VectorType centrality( boost::num_vertices( g ) );
		typename boost::graph_traits< GraphMatrixType >::vertex_iterator vi, vi_end;

		boost::brandes_betweenness_centrality( g,
				centrality_map( boost::get( boost::vertex_centrality, g ) ).
				edge_centrality_map( boost::get( boost::edge_centrality, g ) ).
				weight_map( boost::get(
						boost::edge_weight, g ) ) );

		if ( !relative )
		{
			typename boost::property_map< GraphMatrixType, boost::vertex_centrality_t >::type b = boost::get( boost::vertex_centrality, g );

			for ( tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
				centrality[*vi] = b[*vi];
		}
		else
		{
			boost::relative_betweenness_centrality( g, boost::get( boost::vertex_centrality, g ) );
			typename boost::property_map< GraphMatrixType, boost::vertex_centrality_t >::type r = boost::get( boost::vertex_centrality, g );

			for ( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; ++vi )
				centrality[*vi] = r[*vi];

			// Central point dominance is correct if and only if betweenness centrality is relative!
			if ( dominance )
			{
				double dominance = boost::central_point_dominance( g, boost::get( boost::vertex_centrality, g ) );
				std::cout << "Dominance," << dominance << std::endl;
			}
		}

		WriteVectorToFile( outputFileName, centrality );
	}

	/**
	 * Write vector to filename.
	 */
	void WriteVectorToFile( const std::string& outputFileName, const VectorType& V )
	{
		std::ofstream out( outputFileName.c_str() );

		if ( out.fail() )
		{
			std::cerr << "*** ERROR ***: Not able to write to: " << outputFileName << "!" << std::endl;
			exit( EXIT_FAILURE );
		}

		for ( unsigned int i = 0; i < V.size(); i++ )
			out << V[i] << std::endl;

		out.close();
	}
};

/**
 * Betweenness centrality.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "betweenness centrality", "Calculate betweenness centrality measure for all vertices" );

	std::string inputFileName;
	std::string outputFileName;

	float lowerThreshold = 0;
	float upperThreshold = 1.0;
	bool isUndirected = true;
	bool relative = true;
	bool dominance = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input image: 2D adjacency matrix" ) ->SetRequired( true );

	p.AddArgument( upperThreshold, "threshold" ) ->AddAlias( "thu" ) ->AddAlias( "t" ) ->SetDescription(
			"Threshold; only include paths with weight in (0, threshold] (default: 1.0)" );

	p.AddArgument( lowerThreshold, "threshold-lower" ) ->AddAlias( "thl" ) ->SetDescription(
			"Lower threshold; only include paths with weight > threshold (default: 0)" );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription(
			"Output image: List with (relative) betweenness centralities per vertex" ) ->SetRequired( true );

	p.AddArgument( isUndirected, "undirected" ) ->AddAlias( "u" ) ->SetDescription(
			"Input is a symmetric matrix representing an undirected graph [default: true]" );

	p.AddArgument( relative, "relative" ) ->AddAlias( "r" ) ->SetDescription( "Relative betweenness centrality algorithm [default: true]" );

	p.AddArgument( dominance, "dominance" ) ->AddAlias( "d" ) ->SetDescription( "Print central point dominance [relative only, default: false]" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	if ( isUndirected )
	{
		BC< boost::undirectedS > bc;
		bc.Run( inputFileName, lowerThreshold, upperThreshold, outputFileName, relative, dominance );
	} else
	{
		BC< boost::directedS > bc;
		bc.Run( inputFileName, lowerThreshold, upperThreshold, outputFileName, relative, dominance );
	}

	return EXIT_SUCCESS;
}

