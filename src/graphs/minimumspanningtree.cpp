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
#include <boost/graph/kruskal_min_spanning_tree.hpp>

/**
 * Minimum spanning tree.
 */
class MST
{
public:

	typedef double PixelType;

	typedef boost::property< boost::vertex_index_t, unsigned int, boost::property< boost::vertex_centrality_t, PixelType > >
			VertexPropertyType;
	typedef boost::property< boost::edge_weight_t, PixelType, boost::property< boost::edge_centrality_t, PixelType > > EdgePropertyType;
	typedef boost::adjacency_matrix< boost::undirectedS, VertexPropertyType, EdgePropertyType > GraphMatrixType;
	typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property< boost::edge_weight_t,
			PixelType > > GraphListType;

	typedef itk::Image< PixelType, 2 > ImageType;
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef vnl_vector< PixelType > VectorType;

public:

	/**
	 * Run.
	 */
	void Run( const std::string& inputFileName,
			PixelType lowerThreshold, PixelType upperThreshold, const std::string& outputFileName )
	{
		WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( outputFileName.c_str() );
		writer->SetInput( GetMinimumSpanningTree( inputFileName, lowerThreshold, upperThreshold  ) );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputFileName << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}

	}

protected:

	/**
	 * Return minimum spanning tree as matrix image.
	 */
	ImageType::Pointer GetMinimumSpanningTree( const std::string& filename, PixelType lowerThreshold, PixelType upperThreshold )
	{

		//Example();

		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();
		reader = 0;

		PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
		ImageType::RegionType region = image->GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		int rows = size[0];
		int cols = size[1];

		DataMatrixType data( rows, cols, buffer );

		GraphMatrixType g( rows );

		for ( int i = 0; i < rows; ++i )
		{
			for ( int j = i + 1; j < cols; ++j )
			{
				if ( i == j )
				{
					continue;
				}

				PixelType weight = data( i, j );

				if ( weight <= upperThreshold && weight > 0 && weight > lowerThreshold )
				{
					// invert weight ( e.g. correlation coefficients ... )
					boost::add_edge( i, j, EdgePropertyType( 1.0 / weight ), g );
				}
			}
		}

		// [ 2 ]: Calculate minimum spanning tree ...

		typedef boost::graph_traits< GraphMatrixType >::edge_descriptor EdgeType;
		std::vector< EdgeType > spanningTree;

		// get spanning tree ...
		kruskal_minimum_spanning_tree( g, std::back_inserter( spanningTree ) );

		GraphMatrixType spanningTreeGraph();

		// [ 3 ]: Write minimum spanning tree to image ...

		ImageType::Pointer output = ImageType::New();
		output->CopyInformation( image );
		output->SetRegions( image->GetLargestPossibleRegion() );
		output->Allocate();
		output->FillBuffer( 0 );

		boost::graph_traits< GraphMatrixType >::edge_iterator eiter, eiter_end;
		for ( boost::tie( eiter, eiter_end ) = boost::edges( g ); eiter != eiter_end; ++eiter )
		{
			if ( std::find( spanningTree.begin(), spanningTree.end(), *eiter ) != spanningTree.end() )
			{
				PixelType w = boost::get( boost::edge_weight, g, *eiter );

				ImageType::IndexType index;
				index[0] = source( *eiter, g );
				index[1] = target( *eiter, g );
				output->SetPixel( index, w );

				index[1] = source( *eiter, g );
				index[0] = target( *eiter, g );
				output->SetPixel( index, w );
			}
		}
		return output;
	}
};

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "minimum spanning tree", "Calculate minimum spanning tree on undirected graph" );

	std::string inputFileName;
	std::string outputFileName;

	float lowerThreshold = 0;
	float upperThreshold = 1.0;


	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input image: 2D undirected adjacency matrix" ) ->SetRequired( true );

	p.AddArgument( upperThreshold, "threshold" ) ->AddAlias( "thu" ) ->AddAlias( "t" ) ->SetDescription(
			"Threshold; only include paths with weight in (0, threshold] (default: 1.0)" );

	p.AddArgument( lowerThreshold, "threshold-lower" ) ->AddAlias( "thl" ) ->SetDescription(
			"Lower threshold; only include paths with weight > threshold (default: 0)" );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription(
			"Output image: minimum spanning tree with weights" ) ->SetRequired( true );


	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return -1;
	}

	MST mst;
	mst.Run( inputFileName, lowerThreshold, upperThreshold, outputFileName );

	return EXIT_SUCCESS;
}

