#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkScalarImageToListAdaptor.h"
#include "itkListSampleToHistogramGenerator.h"

/**
 * @author: W.M. Otte
 * @date: 19-11-2009
 * @aim: Determine image histogram from given input image(s).
 * To determine list histograms instead of image histograms, use histogram.
 */
class ImageHistogram
{

public:

	/**
	 * Run.
	 */
	void run( const std::vector< std::string >& inputs, const std::string& output, unsigned int bins, float minValue, float maxValue,
			bool isCumulative, bool reverse, bool percentage )
	{

		// get dimensions...
		unsigned int dims = getDimensions( inputs[0] );

		if ( dims < 2 || dims > 4 )
		{
			std::cerr << "Number of dimensions should be 2, 3 or 4!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// macro for multiple dimension
		#define switchMacro( pixel, dimension ) \
		if ( dims == dimension ) \
		{ \
			process< pixel, dimension >( inputs, output, bins, minValue, maxValue, isCumulative, reverse, percentage ); \
		}

		switchMacro( double, 2 );
		switchMacro( double, 3 );
		switchMacro( double, 4 );

		exit( EXIT_SUCCESS );
	}

private:

	/**
	 * Return cumulative histogram reversed.
	 */
	template< class T >
	void GetCumulativeHistogram( vnl_vector< T >& histogram )
	{
		for ( unsigned int i = 0; i < histogram.size(); i++ )
		{
			vnl_vector< T > tmp = histogram.extract( histogram.size() - i, i );
			histogram[i] = tmp.sum();
		}
	}

	/**
	 * Reverse histogram.
	 */
	template< class T >
	void GetReverseHistogram( vnl_vector< T >& histogram )
	{
		histogram.flip();
	}

	/**
	 * Check dimensions of inputfile. In case of error, Application is terminated.
	 */
	unsigned int getDimensions( const std::string& inputFileName )
	{
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io )
		{
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else
		{
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Process input image.
	 */
	template< class TPixel, unsigned int VDimension >
	void process( const std::vector< std::string >& inputs, const std::string& output, unsigned int bins, float minValue, float maxValue,
			bool cumulative, bool reverse, bool percentage )
	{
		typedef itk::Image< TPixel, VDimension > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;

		typedef itk::ImageFileWriter< ImageType > WriterType;
		typedef itk::Statistics::ScalarImageToListAdaptor< ImageType > AdaptorType;
		typedef TPixel HistogramMeasurementType;
		typedef itk::Statistics::ListSampleToHistogramGenerator< AdaptorType, HistogramMeasurementType > GeneratorType;

		typedef vnl_matrix< TPixel > MatrixType;
		typedef vnl_vector< TPixel > VectorType;

		MatrixType container( bins, inputs.size(), 0.0 );

		std::cout << "Rows (bins): " << container.rows() << std::endl;
		std::cout << "Cols (images): " << container.columns() << std::endl;

		for ( unsigned int i = 0; i < inputs.size(); i++ )
		{
			// read...
			typename ReaderType::Pointer reader = ReaderType::New();
			reader -> SetFileName( inputs[i] );
			reader -> Update();

			// adaptor...
			typename AdaptorType::Pointer adaptor = AdaptorType::New();
			adaptor->SetImage( reader->GetOutput() );

			// generator...
			typename GeneratorType::Pointer generator = GeneratorType::New();
			typename GeneratorType::HistogramType::SizeType size;

			size.Fill( bins );

			generator->SetListSample( adaptor );
			generator->SetNumberOfBins( size );
			generator->SetMarginalScale( 10.0 );

			typename GeneratorType::HistogramType::MeasurementVectorType min;
			typename GeneratorType::HistogramType::MeasurementVectorType max;

			min.Fill( minValue );
			max.Fill( maxValue );

			generator->SetHistogramMin( min );
			generator->SetHistogramMax( max );
			generator->Update();

			typename GeneratorType::HistogramType::ConstPointer histogram = generator->GetOutput();

			VectorType histoVector( bins );

			// update (accumulate) histogram vector...
			for ( unsigned int j = 0; j < bins; j++ )
			{
				histoVector[j] = histogram->GetFrequency( j, 0 );
			}

			if ( cumulative )
			{
				GetCumulativeHistogram( histoVector );
			}
			if ( reverse )
			{
				GetReverseHistogram( histoVector );
			}

			// add to container...
			container.set_column( i, histoVector );
		}

		if ( percentage )
		{
			toPercentage( container );
		}

		writeMatrix( output, "_all.txt", container );
		writeVector( output, "_mean.txt", average( container ) );
		writeVector( output, "_sd.txt", sd( container ) );
	}

	/**
	 * Write vector to filename.
	 */
	template< class TPixel >
	void writeVector( const std::string& outputRoot, const std::string& extra, const vnl_vector< TPixel >& vector )
	{
		std::stringstream ss;
		ss << outputRoot << extra;

		std::ofstream filestream;
		filestream.open( ss.str().c_str() );
		if ( filestream.fail() )
		{
			std::cerr << "ERROR: Not able to write to: " << ss.str() << std::endl;
			exit( EXIT_FAILURE );
		}

		for ( unsigned int i = 0; i < vector.size(); i++ )
		{
			filestream << vector[i] << std::endl;
		}

		filestream.close();
	}

	/**
	 * Write matrix to filename.
	 */
	template< class TPixel >
	void writeMatrix( const std::string& outputRoot, const std::string& extra, const vnl_matrix< TPixel >& matrix )
	{
		std::stringstream ss;
		ss << outputRoot << extra;

		std::ofstream filestream;
		filestream.open( ss.str().c_str() );
		if ( filestream.fail() )
		{
			std::cerr << "ERROR: Not able to write to: " << ss.str() << std::endl;
			exit( EXIT_FAILURE );
		}

		filestream << matrix;;
		filestream.close();
	}

	/**
	 * Average matrix over rows.
	 */
	template< class TPixel >
	vnl_vector< TPixel > average( const vnl_matrix< TPixel >& matrix )
	{
		unsigned int rows = matrix.rows();

		vnl_vector< TPixel > averages( rows, 0.0 );

		for ( unsigned int i = 0; i < rows; i++ )
		{
			vnl_vector< TPixel > row = matrix.get_row( i );
			averages[i] = row.mean();
		}
		return averages;
	}

	/**
	 * SD rows.
	 */
	template< class TPixel >
	vnl_vector< TPixel > sd( const vnl_matrix< TPixel >& matrix )
	{
		unsigned int rows = matrix.rows();
		vnl_vector< TPixel > sds( rows, 0.0 );

		for ( unsigned int i = 0; i < rows; i++ )
		{
			vnl_vector< TPixel > row = matrix.get_row( i );
			sds[i] = sd( row );
		}
		return sds;
	}

	/**
	 * Return container converted to percentage frequencies.
	 */
	template< class TPixel >
	void toPercentage( vnl_matrix< TPixel >& matrix )
	{
		for ( unsigned int i = 0; i < matrix.cols(); i++ )
		{
			vnl_vector< TPixel > column = matrix.get_column( i );

			TPixel sum = column.sum() / 100.0;

			for ( unsigned j = 0; j < column.size(); j++ )
			{
				column[j] /= sum;
			}

			matrix.set_column( i, column );
		}
	}

	/**
	 * Return standard deviation.
	 */
	template< class TPixel >
	TPixel sd( const vnl_vector< TPixel >& vector )
	{
		TPixel size = vector.size();

		if ( size < 2 )
		{
			return 0;
		}

		TPixel mean = vector.mean();

		TPixel sd = 0;
		for ( unsigned int i = 0; i < size; i++ )
		{
			sd += pow( vector[i] - mean, 2 );
		}

		return vcl_sqrt( sd / size );
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::vector< std::string > inputs;
	std::string output;
	int bins = 250;
	float minValue = 0;
	float maxValue = 1;
	bool cumulative = false;
	bool reverse = false;
	bool percentage = true;

	tkd::CmdParser parser( argv[0], "Calculate histogram of one or more input images." );

	parser.AddArgument( inputs, "inputs" ) -> AddAlias( "i" ) -> SetInput( "<strings>" ) -> SetDescription( "Input Image(s) (2D, 3D or 4D)" ) -> SetRequired(
			true );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output Root" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( bins, "bins" ) -> AddAlias( "b" ) -> SetInput( "<int>" ) -> SetDescription( "Number of bins (default: 250)" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( minValue, "min" ) -> AddAlias( "l" ) -> SetInput( "<float>" ) -> SetDescription( "Min value (default: 0)" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( maxValue, "max" ) -> AddAlias( "u" ) -> SetInput( "<float>" ) -> SetDescription( "Max value (default: 1)" ) -> SetRequired(
			false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( cumulative, "cumulative" ) -> AddAlias( "c" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Output cumulative histogram (default: false)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( reverse, "reverse" ) -> AddAlias( "r" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Output reversed (cumulative) histogram (default: false)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( percentage, "percentage" ) -> AddAlias( "p" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Output percentage of total frequencies (default: true)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ImageHistogram imageHistogram = ImageHistogram();

	imageHistogram.run( inputs, output, (unsigned) bins, minValue, maxValue, cumulative, reverse, percentage );
}
