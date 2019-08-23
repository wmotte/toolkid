#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"

#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "tkdCmdParser.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/graph_utility.hpp>

/**
 * Date: 16-11-2009
 */
class CreateSmallWorld {

public:

	typedef double PixelType;
	typedef itk::Image< PixelType, 2 > ImageType;
	typedef vnl_matrix_ref< PixelType > DataMatrixType;
	typedef itk::ImageFileWriter< ImageType > WriterType;

	typedef boost::adjacency_list<> GraphType;
	typedef boost::small_world_iterator< boost::minstd_rand, GraphType > SWGenerator;


	/**
	 * Run.
	 */
	void run( const std::string& output, unsigned int size, unsigned int k, PixelType probability, bool verbose )
	{
		WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( output.c_str() );
		writer->SetInput( GetMatrix( size, k, probability, verbose ) );

		try {
			writer->Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}

	}

protected:

	/**
	 * Return regular matrix.
	 */
	ImageType::Pointer GetMatrix( unsigned int matrixSize, unsigned int k, PixelType probability, bool verbose )
	{
		ImageType::Pointer outputImage = ImageType::New();
		ImageType::SizeType size;
		ImageType::RegionType region;
		size[ 0 ] = matrixSize;
		size[ 1 ] = matrixSize;

		region.SetSize( size );
		outputImage->SetRegions( region );
		outputImage->Allocate();

		boost::minstd_rand generator;

		GraphType g( SWGenerator( generator, matrixSize, k, probability ), SWGenerator(), matrixSize );

		boost::property_map< GraphType, boost::vertex_index_t >::type vertices
			= boost::get( boost::vertex_index, g );

		for( boost::graph_traits< GraphType >::edge_iterator e = edges( g ).first; e != edges( g ).second; ++e )
	    {
			ImageType::IndexType index;
			index[0] = boost::get( vertices, boost::source( *e, g ) );
	    	index[1] = boost::get( vertices, boost::target( *e, g ) );
	    	outputImage->SetPixel( index, static_cast< PixelType >( 1.0 ) );
	    }

		if ( verbose )
			boost::print_graph( g );

		return outputImage;
	}
};

/**
 * Create regular matrix.
 */
int main( int argc, char ** argv ) {
	tkd::CmdParser p( "create_smallworld", "Construct smallworld graph." );

	std::string output;
	int size = 10;
	int k = 2;
	float probability = 0.03;
	bool verbose;

	p.AddArgument( output, "output" ) ->AddAlias( "o" ) ->SetInput( "filename" ) ->SetDescription( "Output matrix file name" ) ->SetRequired(
			true ) -> SetMinMax( 1, 1 );

	p.AddArgument( size, "size" ) ->AddAlias( "s" ) ->SetInput( "uint" ) ->SetDescription( "Number of nodes [ default: 10 ]" );

	p.AddArgument( k, "k" ) ->AddAlias( "k" ) ->SetInput( "uint" ) ->SetDescription( "K [ default: 2 ]" );

	p.AddArgument( probability, "probability" ) ->AddAlias( "p" ) ->SetInput( "float" ) ->SetDescription( "Rewire probability [ default: 0.03 ]" );

	p.AddArgument( verbose, "verbose" ) ->AddAlias( "v" ) ->SetInput( "bool" ) ->SetDescription( "Print graph to screen" );

	if ( !p.Parse( argc, argv ) ) {
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	CreateSmallWorld createSmallWorld;

	createSmallWorld.run( output, size, k, probability, verbose );

	return EXIT_SUCCESS;
}


