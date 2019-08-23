#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkMinimumMaximumImageFilter.h"

#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"

#include "tkdCmdParser.h"

#include "graphCommon.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>

/**
 * Export graph.
 */
class ExportGraph
{
public:
	typedef double PixelType;
	typedef itk::Image< PixelType, 2 > ImageType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef itk::MinimumMaximumImageFilter< ImageType > MinMaxCalcType;

	/**
	 * Run.
	 */
	void Run( const std::string& filename, const std::string& outputFileName,
			bool isPajek, bool isGML, bool isSimilarity, bool isUndirected, const std::string& labelsFileName,
			bool weighted, PixelType threshold, PixelType lineWidth, bool totalDegree, bool averageDegree )
	{
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();

		ImageType::Pointer image = reader->GetOutput();
		reader = 0;

		// alpha of lines is set to weight, so a maximum weight is alpha of 0x00
		// a mimimum weight is alpha of 0xFF.
		// GetImageHexValue return the multiplication value for the weights.
		// (e.g. if the max image value is 3.5, this function returns: 255 / 3.5 ).
		unsigned int hexMultiplicationValue = GetHexImageValue( image );

		PixelType* buffer = image->GetPixelContainer()->GetBufferPointer();
		ImageType::RegionType region = image->GetLargestPossibleRegion();
		ImageType::SizeType size = region.GetSize();
		unsigned int rows = size[0];
		unsigned int cols = size[1];

		DataMatrixType data( rows, cols, buffer );

		if ( isPajek )
			WritePajekGraph( outputFileName + ".net", rows, cols, isUndirected, data );

		if ( isGML )
			WriteGML( outputFileName + ".gml", data, isUndirected, labelsFileName, weighted, hexMultiplicationValue, threshold, "matrix: "
					+ filename, lineWidth, totalDegree, averageDegree );
		
		if ( isSimilarity )
			WriteSimilarity( outputFileName + ".sim", rows, cols, isUndirected, data, labelsFileName );
	}

protected:

	/**
	 * Write Similarity output.
	 */
	void WriteSimilarity( const std::string& outputFileName, 
			unsigned int rows, unsigned int cols, bool isUndirected,
			const DataMatrixType& data, const std::string& labelsFileName )
	{
		// open stream
		std::ofstream out( outputFileName.c_str(), std::ios_base::binary | std::ios_base::out );
		if ( out.fail() )
			std::cerr << "*** ERROR ***: could not write to: " << outputFileName << "!" << std::endl;

		std::vector< std::string > labels;
		bool useLabels = false;
		
		// read labels
		if ( !labelsFileName.empty() )
		{
			graph::Graph< PixelType >::MapType labelMap;
			graph::Graph< PixelType >::GetLabelNames( labelMap, labelsFileName, "\t" );

			graph::Graph< PixelType >::MapType::iterator it = labelMap.begin();
			for( ; it != labelMap.end(); ++it )
			{
				labels.push_back( it->second );
			}
			useLabels = true;
		}	
				
		for ( unsigned int i = 0; i < rows; i++ )
		{
			int start = isUndirected ? i : 0;
			for ( unsigned int j = start; j < cols; ++j )
			{
				PixelType weight = data( i, j );
				
				
				if ( useLabels )
				{
					if ( i == j )
					{
						out << labels.at( i ) << " " << labels.at( j ) << " " << 1 << std::endl; // self-similarity is 1.
					}				
					else if( weight != 0 )
					{
						out << labels.at( i ) << " " << labels.at( j ) << " " << weight << std::endl;
					}
				} 
				else // no labels, so use indices as name ...
				{
					if ( i == j )
					{
						out << i << " " << j << " " << 1 << std::endl; // self-similarity is 1.
					} 
					else if( weight != 0 )
					{
						out << i << " "<< j << " " << weight <<std::endl;
					}
				}
			}
		}

		out.close();
	}
	
	
	
	/**
	 * Return multiplication value to get all image intensities in the range of 0 - 255.
	 */
	unsigned int GetHexImageValue( const ImageType::Pointer image )
	{
		MinMaxCalcType::Pointer filter = MinMaxCalcType::New();
		filter->SetInput( image );
		filter->Update();
		return 255 / static_cast< PixelType > ( filter->GetMaximum() );
	}	
	
	/**
	 * Write graph as GML output
	 *
	 * (visualize with Cytoscape or yEd, etc.
	 * See http://en.wikipedia.org/wiki/Graph_Modelling_Language).
	 *
	 */
	void WriteGML( const std::string& outputFileName, const DataMatrixType& data, bool isUndirected,
			const std::string& labelsFileName, bool weighted,
			unsigned int maxHexValue, PixelType threshold, const std::string& comment,
			PixelType lineWidth, bool totalDegree, bool averageDegree )
	{
		graph::Graph< PixelType >::VectorType nodeLabels;
		std::vector< PixelType > nodeDiameters;
		std::vector< std::string > labels;

		if ( labelsFileName.empty() )
		{
			graph::Graph< PixelType >::Degree( nodeLabels, data, threshold, totalDegree, weighted, false, "", "" );
		}
		else
		{
			graph::Graph< PixelType >::MapType labelMap;
			graph::Graph< PixelType >::GetLabelNames( labelMap, labelsFileName, "\t" );

			graph::Graph< PixelType >::MapType::iterator it = labelMap.begin();
			for( ; it != labelMap.end(); ++it )
			{
				labels.push_back( it->second );
			}
		}

		if( averageDegree )
		{
			std::transform( nodeLabels.begin(), nodeLabels.end(), nodeLabels.begin(),
					std::bind2nd( std::divides< PixelType >(), static_cast< PixelType >( 2.0 ) ) );
		}

		nodeDiameters = NormalizeVector( nodeLabels, 100 );

		// open stream
		std::ofstream out( outputFileName.c_str() );
		if ( out.fail() )
			std::cerr << "*** ERROR ***: could not write to: " << outputFileName << "!" << std::endl;

		WriteGMLRootOpen( out, isUndirected, comment );

		WriteGMLNodes( out, data, nodeDiameters, nodeLabels, labels );

		WriteGMLEdges( out, data, weighted, isUndirected, threshold, maxHexValue, lineWidth );

		WriteGMLRootClose( out );

		out.close();
	}

	/**
	 * Normalize vector so max value is equal to given max.
	 */
	std::vector< PixelType > NormalizeVector( const graph::Graph< PixelType >::VectorType& V, PixelType max )
	{
		std::vector< PixelType > out( V.size() );

		for ( unsigned int i = 0; i < V.size(); i++ )
		{
			out[i] = V[i];
		}

		PixelType maxValue = *( std::max_element( out.begin(), out.end() ) );
		PixelType multiply = static_cast< PixelType > ( max ) / maxValue;

		std::transform( out.begin(), out.end(), out.begin(), std::bind2nd( std::multiplies< PixelType >(), multiply ) );

		return out;
	}

	/**
	 * Write root tags in GML output.
	 *
	 */
	void WriteGMLRootOpen( std::ofstream& out, bool isUndirected, const std::string& label )
	{
		std::string indent = "\t";

		out << "Creator " << "\"wim@invivonmr.uu.nl (InvivoNMR UMC Utrecht)\"" << std::endl;
		out << "Version \"2.7\"" << std::endl;
		out << "graph [" << std::endl;
		out << indent << "hierarchic " << 1 << std::endl;
		out << indent << "label \"" << label << "\"" << std::endl;
		out << indent << "directed " << ( isUndirected ? 0 : 1 ) << std::endl;
	}

	/**
	 * Close graph opject.
	 */
	void WriteGMLRootClose( std::ofstream& out )
	{
		out << "]" << std::endl;
	}

	/**
	 * Write GML Node.
	 */
	void WriteGMLNode( std::ofstream& out, int id, PixelType x, PixelType y,
			PixelType diameter, PixelType degree, const std::string& label )
	{
		if ( diameter == 0 )
		{
			std::cerr << "*** ERROR ***: writing GML node. Node diameter should be non-zero!" << std::endl;
			out.close();
			exit( EXIT_FAILURE );
		}

		PixelType dropShadowOffset = diameter / static_cast< PixelType > ( 12.0 );
		unsigned int fontsize = ( diameter / static_cast< PixelType > ( 2 ) ) + 1;

		std::string indent = "\t";
		std::string nindent = "\t\t";
		std::string gindent = "\t\t\t";

		out << indent << "node [" << std::endl;
		out << nindent << "id " << id << std::endl;
		out << nindent << "graphics [" << std::endl;

		out << gindent << "x " << x << std::endl;
		out << gindent << "y " << y << std::endl;

		out << gindent << "w " << diameter << std::endl;
		out << gindent << "h " << diameter << std::endl;

		out << gindent << "fill " << "\"#FFFF00\"" << std::endl;
		out << gindent << "fill2 " << "\"#FF0000\"" << std::endl;
		out << gindent << "type " << "\"ellipse\"" << std::endl;
		out << gindent << "outline " << "\"#000000\"" << std::endl;

		out << gindent << "dropShadowColor " << "\"#00000040\"" << std::endl;
		out << gindent << "dropShadowOffsetX " << dropShadowOffset << std::endl;
		out << gindent << "dropShadowOffsetY " << dropShadowOffset << std::endl;

		// close graphics
		out << nindent << "]" << std::endl;

		out << nindent << "LabelGraphics [" << std::endl;

		if ( ! label.empty() )
		{
			fontsize = ( diameter / static_cast< PixelType > ( 7 ) ) + 1;
			out << gindent << "fontSize " << fontsize << std::endl;
			out << gindent << "text " << "\"" << label << "\"" << std::endl;
		}
		else
		{
			out << gindent << "fontSize " << fontsize << std::endl;
			std::streamsize old = out.precision(); // save current
			out << gindent << "text " << "\"" << std::setprecision( 4 ) << degree << "\"" << std::endl;
			out.precision( old ); // restore
		}

		out << gindent << "outline \"#808080\"" << std::endl;
		out << gindent << "fill \"#FFFFFF80\"" << std::endl;
		out << gindent << "fontName \"Dialog\"" << std::endl;
		out << gindent << "anchor \"c\"" << std::endl;
		out << nindent << "]" << std::endl;

		// close node
		out << indent << "]" << std::endl;
	}

	/**
	 * Write GML Edge.
	 */
	void WriteGMLEdge( std::ofstream& out, int root_index, int target,
			int source, unsigned int width, const std::string& label, PixelType lineWidth )
	{
		if ( width > 255 )
		{
			std::cerr << "*** WARNING ***: overflow -> setting GML transparence to FF!" << std::endl;
			width = 255;
		}

		std::string indent = "\t";
		std::string eindent = "\t\t";
		std::string gindent = "\t\t\t";

		out << indent << "edge [" << std::endl;
		out << eindent << "source " << source << std::endl;
		out << eindent << "target " << target << std::endl;

		out << eindent << "graphics [" << std::endl;

		// argb -> convert weight to transparence ...
		out << gindent << "width " << lineWidth << std::endl;
		out << gindent << "fill " << "\"#000000" << std::hex << width << "\"" << std::dec << std::endl;

		out << gindent << "type " << "\"line\"" << std::endl;

		out << gindent << "Line [ " << std::endl;
		out << gindent << "]" << std::endl;

		// close graphics
		out << eindent << "]" << std::endl;
		out << eindent << "label \"" << label << "\"" << std::endl;

		// close edge
		out << indent << "]" << std::endl;
	}

	/**
	 * Write GML Nodes.
	 *
	 * Return last root_index.
	 */
	void WriteGMLNodes( std::ofstream& out, const DataMatrixType& data,
			const std::vector< PixelType >& diameters,
			const graph::Graph< PixelType >::VectorType& degrees,
			const std::vector< std::string > labels )
	{
		unsigned int rows = data.rows();
		bool useDiameter = false;

		if ( !diameters.empty() )
		{
			if ( diameters.size() != rows )
			{
				std::cerr << "*** ERROR ***: "
					"Could not write GML nodes: number of diameters /= matrix rows!" << std::endl;
				out.close();
				exit( EXIT_FAILURE );
			}
			useDiameter = true;
		}

		// x,y,z coordinates set to id, for the moment ...
		for ( unsigned int i = 0; i < rows; i++ )
		{
			std::stringstream ss;
			ss << i + 1;

			if ( ! labels.empty() )
			{
				std::string label = labels[i];
				WriteGMLNode( out, i + 1, i + 1, i + 1, 100.0, 1.0, label );
			}
			else if ( useDiameter )
			{
				PixelType diameter = diameters[i];
				PixelType degree = degrees[i];

				// if diameter == 0, it is an isolated node -> remove!
				if ( diameter != 0 )
				{
					WriteGMLNode( out, i + 1, i + 1, i + 1, diameter, degree, ss.str() );
				}
			}
			else
			{
				WriteGMLNode( out, i + 1, i + 1, i + 1, 1.0, 1.0, ss.str() );
			}
		}
	}

	/**
	 * Write GML Edges.
	 *
	 * Start root_indices with given root_index start value.
	 * Write edge values if weighted is true.
	 */
	void WriteGMLEdges( std::ofstream& out, const DataMatrixType& data, bool weighted, bool isUndirected, PixelType threshold,
			unsigned int maxHexValue, PixelType lineWidth )
	{
		unsigned int rows = data.rows();
		unsigned int cols = data.columns();

		int root_index = 0;

		for ( unsigned int i = 0; i < rows; i++ )
		{
			int start = isUndirected ? i + 1 : 0;
			for ( unsigned int j = start; j < cols; ++j )
			{
				if ( i == j )
				{
					continue;
				}

				PixelType weight = data( i, j );

				if ( weight > threshold )
				{
					if ( weighted )
						WriteGMLEdge( out, root_index + 1, j + 1, i + 1, weight * maxHexValue, "", lineWidth );
					else
						WriteGMLEdge( out, root_index + 1, j + 1, i + 1, 255, "", lineWidth ); // all lines black
					root_index++;
				}
			}
		}
	}

	/**
	 * Write Pajek graph.
	 */
	void WritePajekGraph( const std::string& outputFileName, unsigned int rows, unsigned int cols, bool isUndirected,
			const DataMatrixType& data )
	{
		// open stream
		// open stream
		std::ofstream out( outputFileName.c_str(), std::ios_base::binary | std::ios_base::out );
		if ( out.fail() )
			std::cerr << "*** ERROR ***: could not write to: " << outputFileName << "!" << std::endl;

		// write Pajek Vertices
		out << "*Vertices " << rows << "\r\n";

		for ( unsigned int i = 0; i < rows; i++ )
		{
			std::stringstream ss;
			ss << "node_" << ( i + 1 );
			unsigned int pi = i + 1;
			WritePajekNode( out, pi, ss.str() );
		}

		// write Pajek Edges
		out << "*Edges" << "\r\n";

		for ( unsigned int i = 0; i < rows; i++ )
		{
			int start = isUndirected ? i + 1 : 0;
			for ( unsigned int j = start; j < cols; ++j )
			{
				if ( i == j )
				{
					continue;
				}

				PixelType weight = data( i, j );

				if ( weight != 0 )
				{
					unsigned int pi = i + 1;
					unsigned int pj = j + 1;
					WritePajekEdge( out, pi, pj, weight );
				}

			}
		}

		out.close();
	}

	/**
	 * Write node to output stream.
	 */
	void WritePajekNode( std::ofstream& out, unsigned int i, const std::string& name )
	{
		out << i << " " << name << "\r\n";
	}

	/**
	 * Write edge to output stream.
	 */
	void WritePajekEdge( std::ofstream& out, unsigned int i, unsigned int j, PixelType w )
	{
		out << i << " " << j << " " << w << "\r\n";
	}

};

/**
 * Export Matrix -> Pajek, Graphviz, GML.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( "export", "Export graph to Pajek, GML or Similarity output" );

	std::string inputFileName;
	std::string outputFileName;

	std::string labelsFileName;

	bool isPajek = false;
	bool isSimilarity = false;
	bool isGML = false;
	
	bool isUndirected = true;
	bool weighted = true;

	float threshold = 0.0;
	float lineWidth = 2.0;

	bool totalDegree = false;

	bool averageDegree = false;

	p.AddArgument( inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input image: 2D adjacency matrix" ) ->SetRequired( true );

	p.AddArgument( outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output root name" ) ->SetRequired( true );

	p.AddArgument( isPajek, "pajek" ) ->AddAlias( "p" ) ->SetDescription( "Export to Pajek [default: false]" );

	p.AddArgument( isSimilarity, "sim" ) ->AddAlias( "s" ) ->SetDescription( "Export to Similarity [default: false]" );
	
	p.AddArgument( isGML, "gml" ) ->AddAlias( "m" ) ->SetDescription( "Export to GML [default: false]" );

	p.AddArgument( totalDegree, "total-degree" ) ->AddAlias( "td" ) ->SetDescription( "Total degree instead of out-degree [default: false]" );

	p.AddArgument( averageDegree, "average-degree" ) ->AddAlias( "ad" ) ->SetDescription( "Divide total-degree with factor 2 to get average in-degree and out-degree [default: false]" );

	p.AddArgument( isUndirected, "undirected" ) ->AddAlias( "u" ) ->SetDescription(
			"Input is a symmetric matrix representing an undirected graph [default: true]" );

	p.AddArgument( labelsFileName, "labels" ) ->AddAlias( "d" ) ->SetDescription( "Insert text file labels in nodes [GML only; if no file is specified, node degree is added as label]" );

	p.AddArgument( weighted, "weighted" ) ->AddAlias( "w" ) ->SetDescription( "Use weighted edges in GML [default: true]" );

	p.AddArgument( threshold, "threshold" ) ->AddAlias( "t" ) ->SetDescription( "Exclude edges lower or equal to threshold [default: 0.0]" );

	p.AddArgument( lineWidth, "lineWidth" ) ->AddAlias( "l" ) ->SetDescription( "Edges line width [default: 2.0]" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ExportGraph exp;
	exp.Run( inputFileName, outputFileName, isPajek, isGML, isSimilarity, isUndirected, labelsFileName, weighted, threshold, lineWidth, totalDegree, averageDegree );

	return EXIT_SUCCESS;
}
