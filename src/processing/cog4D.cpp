#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPoint.h"

#include "vnl/algo/vnl_fft_1d.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

/**
 * Wim Otte wim@invivonmr.uu.nl
 *
 * Date: 10-11-2009
 *
 * In this class the center of gravity is calculated for each 3D volume of 4D input;
 * 1: Average 3D point calculated
 * 2: Euclidean distance for all cogs to average 3D point calculated
 * 3: Normalized using standard-deviation of euclidean distance and euclidean mean
 * 4: Outliers determined using given threshold.
 */

class Cof4D
{

public:

	typedef float PixelType;
	typedef unsigned int BinaryPixelType;

	/**
	 * Run cof4D application.
	 */
	void run( const std::string& input, const std::string& output, double threshold,
				bool useNormalized, bool useNeighborhood, unsigned int neighborhoodSize,
					PixelType highPass, PixelType lowPass )
	{
		// get dimensions...
		unsigned int dims = getDimensions( input );

		if ( dims == 4 )
		{
			process( input, output, (PixelType) threshold, useNormalized, useNeighborhood, neighborhoodSize, highPass, lowPass );
		}
		else
		{
			std::cerr << "Number of input dimensions should be 4!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

private:
	typedef itk::Point< PixelType, 3 > Point3DType;
	typedef itk::Image< PixelType, 3 > Image3DType;
	typedef itk::Image< PixelType, 4 > Image4DType;
	typedef itk::ImageFileReader< Image4DType > Reader4DType;
	typedef itk::ImageFileWriter< Image4DType > Writer4DType;
	typedef itk::ExtractImageFilter< Image4DType, Image3DType > ExtractImageFilter3DType;
	typedef itk::ImageRegionIteratorWithIndex< Image3DType > Iterator3DType;
	typedef vnl_matrix_ref< PixelType > ImageMatrixType;
	typedef vnl_vector< PixelType > VectorType;
	typedef vnl_vector< BinaryPixelType > BinaryVectorType;

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
		}
		else
		{
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Process 4D input.
	 */
	void process( const std::string& input, const std::string& output,
			PixelType threshold, bool useNormalized, bool useNeighborhood,
				unsigned int neighborhoodSize, PixelType highPass, PixelType lowPass )
	{
		// [ 1 ]: get cogs.
		std::vector< Point3DType > cogs = getCogs( input );

		// [ 2 ]: calculate mean point.
		Point3DType meanPoint = getMeanPoint( cogs );

		// [ 3 ]: calculate Euclidean distances.
		VectorType distances = getDistances( cogs, meanPoint );

		// [ 4 ]: combine bandpass with detrending.
		FourierFilter( distances, getSampleSpacing( input ), highPass, lowPass );
		RemoveTrend( distances, 0, false );

		// [ 5 ]: get mean.
		PixelType mean = getMean( distances );

		// [ 6 ]: get sd.
		PixelType sd = getSD( distances );

		// [ 7 ]: normalize cogs -> ncogs.
		VectorType ncogs = normalizeDistances( distances, mean, sd );

		// [ 8 ]: get normalized mean
		PixelType nMean = getMean( ncogs );

		// [ 9 ]: get normalized sd
		PixelType nSD = getSD( ncogs );

		// [ 10 ]: get outliers based on threshold.
		BinaryVectorType outliers;

		if ( useNormalized )
		{
			outliers = getOutliers( ncogs, threshold, useNeighborhood, neighborhoodSize );
		}
		else
		{
			outliers = getOutliers( distances, threshold, useNeighborhood, neighborhoodSize );
		}

		// [ 11 ]: get percentage outliers.
		PixelType percentageOutliers = getPercentageOutliers( outliers );

		// [ 12 ]: output_parameters.txt
		writeOutputParameters( output, meanPoint, mean, sd, percentageOutliers, threshold, nMean, nSD, useNormalized );

		// [ 13 ]: output_outliers.txt
		writeVector( output, "_outliers.txt", outliers );

		// [ 14 ]: output_distances.txt
		writeVector( output, "_distances.txt", distances );

		// [ 15 ]: output_ncogs.txt
		writeVector( output, "_ncogs.txt", ncogs );

		// [ 16 ]: output.nii.gz
		writeOutput( input, output, outliers );
	}

	/**
	 * Fourier filter vector.
	 */
	VectorType FourierFilter( VectorType distances, PixelType sampleSpacing, PixelType highPass, PixelType lowPass )
	{
		unsigned int numberOfTimePoints = distances.size();
		int bufferLength = pow( 2, vcl_ceil( vcl_log( numberOfTimePoints ) / vcl_log( 2. ) ) );
		float sampleFactor = static_cast< float >( bufferLength ) * sampleSpacing;
		int lowerIndex1 = static_cast< int >( vcl_ceil( highPass * sampleFactor ) );
		int higherIndex1 = static_cast< int >( vcl_floor( lowPass * sampleFactor ) );
		int lowerIndex2 = bufferLength - higherIndex1 - 1;
		int higherIndex2 = bufferLength - lowerIndex1 - 1;

		vnl_vector< vcl_complex< float > > buffer( bufferLength );
		buffer.fill( itk::NumericTraits< vcl_complex< float > >::Zero );

	    float a = 0;
	    float b = 0;
	    float mean = 0;

	    // FITLINE...
	    FitLine( distances, 0, a, b, mean );

	    for( unsigned int j = 0; j < numberOfTimePoints; ++j )
	    {
	    	float onLine = a * static_cast< float >( j ) + b;
	    	buffer( j ) = vcl_complex< float >( distances( j ) - onLine, 0 );
	    }

	    // zero-filling to power of 2
	    for( int j = numberOfTimePoints; j < bufferLength; ++j )
	    {
	    	buffer( j ) = 0;
	    }

	    // forward transform
	    Fourier< float >( reinterpret_cast< float* >( buffer.data_block() ), buffer.size(), 1 );

	    for( int j = 0; j < bufferLength; ++j )
	    {
	    	if ( j >= lowerIndex1 && j <= higherIndex1 )
	        {
	    		continue;
	        }

	    	if ( j >= lowerIndex2 && j <= higherIndex2 )
	        {
	    		continue;
	        }

	    	buffer( j ) = 0;
	    }

	    // backward transform
	    Fourier< float >( reinterpret_cast< float* >( buffer.data_block() ), buffer.size(), -1 );

	    for( unsigned int j = 0; j < numberOfTimePoints; ++j )
	    {
	    	const vcl_complex< float >& p = buffer( j );
	    	distances( j ) = p.real();
	    }
		return distances;
	}

	/**
	 * Fit linear line.
	 * (c) Press et al., Numerical Recipes 3rd ed., p.784 (Overgenomen van Kajo).
	 */
	void FitLine( const VectorType& data, int ignore, PixelType& a, PixelType& b, PixelType& mean )
	{
		int numberOfPoints = data.size();

		PixelType sx = 0;
	    PixelType sy = 0;
	    PixelType st2 = 0;

	    for( int i = ignore; i < ( numberOfPoints - ignore ); ++i )
	    {
	    	sx += static_cast< PixelType >( i );
	    	sy += data( i );
	    }

	    PixelType sxoss = sx / static_cast< PixelType >( numberOfPoints - 2 * ignore);

	    for( int i = ignore; i < ( numberOfPoints - ignore ); ++i )
	    {
	    	PixelType t = static_cast< PixelType >( i ) - sxoss;
	    	st2 += t * t;
	    	a += t * data( i );
	    }
	    a /= st2;
	    b = ( sy - sx * a ) / static_cast< PixelType >( numberOfPoints - 2 * ignore );
	    mean = sy / static_cast< PixelType >( numberOfPoints - 2 * ignore );
	}

	/**
	 * Swap data.
	 */
	template< class T >
	void Swap( T& a, T& b )
	{
		T tmp = a;
		a = b;
		b = tmp;
	}

	/**
	 * Fourier transform.
	 * (c) Press et al., Numerical Recipes 3rd ed., p.612 (overgenomen van Kajo)
	 */
	template< class T >
	void Fourier( T* data, const int& n, const int& isign )
	{
		int nn, mmax, m, j, istep, i;
		T wtemp, wr, wpr, wpi, wi, theta, tempr, tempi, norm;

		if ( n < 2 || n & ( n - 1 ) )
	    {
			throw( "n must be a power of 2" );
	    }

		nn = n << 1;
		j = 1;

		for( i = 1; i < nn; i += 2 )
	    {
			if ( j > i )
			{
				Swap< T >( data[ j - 1 ], data[ i - 1 ] );
				Swap< T >( data[ j ], data[ i ] );
			}
			m = n;
			while( m >= 2 && j > m )
			{
				j -= m;
				m >>= 1;
			}
			j += m;
	    }

		mmax = 2;

		while( nn > mmax )
	    {
			istep = mmax << 1;
			theta = static_cast< T >( isign ) * ( 6.28318530717959 / mmax );
			wtemp = sin( 0.5 * theta );
			wpr = -2.0 * wtemp * wtemp;
			wpi = sin( theta );
			wr = 1.0;
			wi = 0.0;

			for( m = 1; m < mmax; m +=2 )
			{
				for( i = m; i <= nn; i += istep )
				{
					j = i + mmax;
					tempr = wr * data[ j - 1 ] - wi * data[ j ];
					tempi = wr * data[ j ] + wi * data[ j - 1 ];
					data[ j - 1 ] = data[ i - 1 ] - tempr;
					data[ j ] = data[ i ] - tempi;
					data[ i - 1 ] += tempr;
					data[ i ] += tempi;
				}
				wtemp = wr;
				wr = wtemp * wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
			}
			mmax = istep;
	    }

		if ( isign == -1 )
	    {
			norm = 1. / n;
			for( i = 0; i < nn; ++i )
			{
				data[ i ] *= norm;
			}
	    }
	}

	/**
	 * Remove linear trend.
	 */
	void RemoveTrend( VectorType& signal, int ignore = 0, bool addMean = false )
	{
		PixelType a = 0;
		PixelType b = 0;
		PixelType mean = 0;
	    this->FitLine( signal, ignore, a, b, mean );

	    const int n = signal.size();
	    for( int i = 0; i < n; ++i )
		{
	    	signal( i ) = signal( i ) - a * static_cast< PixelType >( i ) - b;
	    	if( addMean )
	        {
	    		signal( i ) += mean;
	        }
	    }
	}

	/**
	 * Write output 4D image with outlier 3D volumes removed.
	 */
	template<class T>
	void writeOutput( const std::string& input, const std::string& output, vnl_vector<T>& outliers )
	{
		// get output file name...
		std::stringstream ss;
		ss << output << ".nii.gz";
		std::string outputFileName = ss.str();

		// read input volume...
		Reader4DType::Pointer reader = Reader4DType::New();
		reader -> SetFileName( input );
		Image4DType::Pointer inputImage = reader->GetOutput();
		reader -> Update();

		// 4D region...
		Image4DType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
		Image4DType::IndexType index = region.GetIndex();
		Image4DType::SizeType size = region.GetSize();

		// initialize new 4D image...
		Image4DType::Pointer outputImage = Image4DType::New();
		outputImage -> CopyInformation( inputImage );
		Image4DType::RegionType outputRegion = region;

		// get nr. of output/exclusion volumes.
		unsigned int totalVolumes = size[3];
		unsigned int totalExclusionVolumes = outliers.sum();

		std::cout << "total exclusion volumes: " << totalExclusionVolumes << std::endl;

		size[3] = size[3] - totalExclusionVolumes;
		outputRegion.SetIndex( index );
		outputRegion.SetSize( size );
		outputImage -> SetRegions( outputRegion );
		outputImage -> Allocate();

		// for each volume...
		unsigned int j = 0;

		// prepare for iteration...
		Image4DType::IndexType inputIndex = region.GetIndex();
		Image4DType::SizeType inputSize = region.GetSize();
		inputSize[3] = 1;

		Image4DType::IndexType outputIndex = outputRegion.GetIndex();
		Image4DType::SizeType outputSize = outputRegion.GetSize();
		outputSize[3] = 1;

		for( unsigned int i = 0; i < totalVolumes; i++ )
		{
			if ( outliers[i] != 1 )
			{
				Image4DType::RegionType inputRegion;
				Image4DType::RegionType outputRegion;

				inputIndex[3] = i;
				outputIndex[3] = j;

				inputRegion.SetIndex( inputIndex );
				inputRegion.SetSize( inputSize );

				outputRegion.SetIndex( outputIndex );
				outputRegion.SetSize( outputSize );

				itk::ImageRegionConstIteratorWithIndex< Image4DType > inputIterator( inputImage, inputRegion );
				itk::ImageRegionIteratorWithIndex< Image4DType > outputIterator( outputImage, outputRegion );

				for ( ; !inputIterator.IsAtEnd(), !outputIterator.IsAtEnd(); ++inputIterator, ++outputIterator )
				{
					Image4DType::IndexType pixelIndex = inputIterator.GetIndex();
					outputIterator.Set( inputImage -> GetPixel( pixelIndex ) );
				}

				// [ 3 ] update output 4D index...
				j++;
			}
		}

		// Write output image.
		Writer4DType::Pointer writer = Writer4DType::New();
		writer -> SetFileName( outputFileName.c_str() );
		writer -> SetInput( outputImage );

		try
		{
			writer -> Update();
		}
		catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputFileName << std::endl;
		}
	}


	/**
	 * Write vector to text.
	 */
	template <class T>
	void writeVector( const std::string& output,
											const std::string& extra,
												const vnl_vector<T>& vector )
	{
		std::stringstream ss;
		ss << output << extra;

		std::ofstream filestream;
		filestream.open( ss.str().c_str() );
		if ( filestream.fail() ) {
			std::cerr << "ERROR: Not able to write to: " << ss.str() << std::endl;
			exit( EXIT_FAILURE );
		}

		// for each potential outliers...
		for( unsigned int i = 0; i < vector.size(); i++ )
		{
			filestream << vector[i] << std::endl;
		}

		filestream.close();
	}

	/**
	 * Write parameters.
	 */
	void writeOutputParameters( const std::string& output, Point3DType meanPoint,
									PixelType mean, PixelType sd, PixelType percentageOutliers,
										PixelType threshold, PixelType nMean, PixelType nSD, bool useNormalized )
	{
		std::stringstream ss;
		ss << output << "_parameters.txt";

		std::ofstream filestream;
		filestream.open( ss.str().c_str() );
		if ( filestream.fail() ) {
			std::cerr << "ERROR: Not able to write to: " << ss.str() << std::endl;
			exit( EXIT_FAILURE );
		}

		filestream << "Distance_Mean, " << mean << std::endl;
		filestream << "Distance_SD, " << sd << std::endl;
		filestream << "Distance_Threshold, " << threshold << std::endl;
		filestream << "Distance_Percentage_Outliers, " << percentageOutliers << std::endl;
		filestream << "Distance_Normalized_Mean, " << nMean << std::endl;
		filestream << "Distance_Normalized_SD, " << nSD << std::endl;
		filestream << "Average_COG, " << meanPoint << std::endl;
		filestream << "Used_Normalized_Data, " << useNormalized << std::endl;

		filestream.close();
	}

	/**
	 * Get percentage outliers.
	 */
	template<class T>
	PixelType getPercentageOutliers( const vnl_vector<T>& outliers )
	{
		return ( static_cast< PixelType >( outliers.sum() ) / static_cast< PixelType>( outliers.size() ) ) * 100.0;
	}

	/**
	 * Return integer vector with 1 for outlier and 0 for non-outlier volumes.
	 */
	template<class T>
	BinaryVectorType getOutliers( const vnl_vector<T>& distances, PixelType threshold, bool useNeighborhood, unsigned int neighborhoodSize )
	{
		BinaryVectorType outliers( distances.size() );

		// for each normalized distance...
		unsigned int size = distances.size();

		for ( unsigned int i = 0; i < size; i++ )
		{
			if ( fabs(distances[i]) > threshold )
			{
				// mark volume as outlier...
				outliers[i] = 1;

				// optionally, mark neighboring voxel(s) as outliers...
				if ( useNeighborhood )
				{
					if ( i > neighborhoodSize - 1 )
					{
						for ( unsigned int j = 1; j <= neighborhoodSize; j++ )
						{
							outliers[i-j] = 1;

						}
					}

					if ( i < size - neighborhoodSize )
					{
						for ( unsigned int j = 1; j <= neighborhoodSize; j++ )
						{
							outliers[i+j] = 1;
						}

						//  update i with its positive neighbors...
						i += neighborhoodSize;
					}
				}
			}
		}
		return outliers;
	}

	/**
	 * Return mean.
	 */
	template<class T>
	PixelType getMean( const vnl_vector<T>& distances )
	{
		return distances.mean();
	}

	/**
	 * Return standard deviation.
	 */
	template<class T>
	PixelType getSD( const vnl_vector<T>& distances )
	{
		PixelType size = distances.size();

		if ( size < 2 )
		{
			return 0;
		}

		PixelType mean = distances.mean();

		PixelType sd = 0;
		for( unsigned int i = 0; i < size; i++ )
		{
			sd += pow( distances[i] - mean, 2 );
		}

        return vcl_sqrt( sd / size );
     }

	/**
	 * Normalize distances using mean and sd transform. N = m - M / sigma.
	 */
	template<class T>
	vnl_vector<T> normalizeDistances( const vnl_vector<T>& distances, PixelType mean, PixelType sd )
	{
		vnl_vector<T> nDistances( distances.size() );

		for ( unsigned int i = 0; i < distances.size(); i++ )
		{
			nDistances[i] = ( distances[i] - mean ) / sd;
		}
		return nDistances;
	}

	/**
	 * Return vector with Euclidean distances.
	 */
	VectorType getDistances( const std::vector< Point3DType >& cogs, Point3DType mean )
	{
		VectorType distances( cogs.size() );

		// for each point...
		for ( unsigned int i = 0; i < cogs.size(); i++ )
		{
			Point3DType cog = cogs[i];
			PixelType distance = 0;

			// for each dimension...
			for ( unsigned int j = 0; j < 3; j++ )
			{
				distance += pow( cog[j] - mean[j], 2 );
			}

			distances[i] = sqrt( distance );
		}
		return distances;
	}

	/**
	 * Return average point.
	 */
	Point3DType getMeanPoint( const std::vector< Point3DType >& cogs )
	{
		Point3DType average;
		average[0] = 0;
		average[1] = 0;
		average[2] = 0;

		for( unsigned int i=0; i < cogs.size(); i++ )
		{
			Point3DType cog = cogs[i];
			average[0] =( average[0] + cog[0] ) / 2.0;
			average[1] =( average[1] + cog[1] ) / 2.0;
			average[2] =( average[2] + cog[2] ) / 2.0;
		}
		return average;
	}

	/**
	 * Return sample spacing (time) of given 4D-volume.
	 */
	PixelType getSampleSpacing( const std::string& input )
	{
		Reader4DType::Pointer reader = Reader4DType::New();
		reader -> SetFileName( input );
		reader -> Update();
		Image4DType::RegionType region = reader -> GetOutput() -> GetLargestPossibleRegion();
		Image4DType::SizeType size = region.GetSize();

		return size[3];
	}

	/**
	 * Return vector with cogs for each 3D volume.
	 */
	std::vector< Point3DType > getCogs( const std::string& input )
	{

		// read...
		Reader4DType::Pointer reader = Reader4DType::New();
		reader -> SetFileName( input );
		reader -> Update();

		// region...
		Image4DType::RegionType region = reader -> GetOutput() -> GetLargestPossibleRegion();
		Image4DType::IndexType index = region.GetIndex();
		Image4DType::SizeType size = region.GetSize();

		unsigned int totalVolumes = size[3];
		size[3] = 0;

		std::vector< Point3DType > cogs;

		// for each volume...
		for( unsigned int i = 0; i < totalVolumes; i++ )
		{
			index[3] = i;
			region.SetSize( size );
			region.SetIndex( index );

			ExtractImageFilter3DType::Pointer filter = ExtractImageFilter3DType::New();

			filter -> SetInput( reader -> GetOutput() );
			filter -> SetExtractionRegion( region );
			Image3DType::Pointer volume = filter -> GetOutput();
			filter -> Update();
			Point3DType cog = getCog( volume );
			cogs.push_back( cog );
		}

		return cogs;
	}

	/**
	 * Return center-of-gravity for 3D volume.
	 */
	Point3DType getCog( const Image3DType::Pointer volume )
	{
		PixelType posx = 0;
		PixelType posy = 0;
		PixelType posz = 0;
		PixelType total_masses = 0;

		Image3DType::RegionType region = volume -> GetLargestPossibleRegion();

		Image3DType::SpacingType spacing = volume -> GetSpacing();

		Iterator3DType it( volume, region );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			Image3DType::IndexType index = it.GetIndex();
			PixelType mass = volume -> GetPixel( index );

			posx += mass * index[0];
			posy += mass * index[1];
			posz += mass * index[2];

			total_masses += mass;
		}

		Point3DType point;
		point[0] = ( posx * spacing[0] ) / total_masses;
		point[1] = ( posy * spacing[1] ) / total_masses;
		point[2] = ( posz * spacing[2] ) / total_masses;

		return point;
	}
};

// ***********************************************************************************

/**
 * Main.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::string input;
	std::string output;
	double threshold = 0.1;
	double highPass = 0.01;
	double lowPass = 1000.0;
	bool useNormalized = false;
	bool useNeighborhood = true;
	int neighborhoodSize = 4;


	tkd::CmdParser parser( argv[0], "Calculate (normalized) center-of-gravity for 3D volumes\n "
									"and discard all volumes (with neighbors) above threshold \n"
									"(n.b. data is detrended first)" );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "4D Input File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output root" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( threshold, "threshold" ) -> AddAlias( "t" ) -> SetInput( "<double>" ) -> SetDescription(
			"(Normalized) center-of-gravity threshold (Euclidean distance, default: 0.1 mm)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( useNormalized, "normalized" ) -> AddAlias( "z" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Use Normalized calculation instead of raw (mm) cog-values (default: false)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( useNeighborhood, "neighborhood" ) -> AddAlias( "n" ) -> SetInput( "<bool>" ) -> SetDescription(
			"Use Neighborhood deletion (mark neighboring volumes as outliers if volume is out of range (default: true)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( neighborhoodSize, "neighborhood-size" ) -> AddAlias( "ns" ) -> SetInput( "<uint>" ) -> SetDescription(
			"Neighborhood size (number of neighboring volumes (default: 4; n.b. total neighbors is 8 (positive and negative))" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( highPass, "high-pass" ) -> AddAlias( "h" ) -> SetInput( "<float>" ) -> SetDescription(
			"High-pass filter (default: 0.01)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	parser.AddArgument( lowPass, "low-pass" ) -> AddAlias( "l" ) -> SetInput( "<float>" ) -> SetDescription(
			"Low-pass filter (default: 1000.0)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Cof4D cof4d = Cof4D();

	cof4d.run( input, output, threshold, useNormalized, useNeighborhood, (unsigned) neighborhoodSize, highPass, lowPass );
}
