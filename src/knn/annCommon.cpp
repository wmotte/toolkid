#include "annCommon.h"

/**
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC Utrecht, NL.
 * @date: 26-12-2009
 *
 * Function implementations of ann common.
 */
namespace ann
{

	/**
	 * Convert vnl matrix to std matrix.
	 */
	template< class ValueType >
	void ann::Ann< ValueType >::ConvertVNL2STD( const vnl_matrix< ValueType >& G, MatrixType& stdG )
	{
		int rows = G.rows();
		int cols = G.cols();

		stdG = MatrixType( rows, VectorType( cols, 0 ) );
		for( int i = 0; i < rows; i++ )
			for( int j = 0; j < cols; j++ )
				stdG[i][j] = G( i, j );
	}

	/**
	 * Convert vnl vector to std vector.
	 */
	template< class ValueType >
	void ann::Ann< ValueType >::ConvertVNL2STD( const vnl_vector< ValueType >& G, VectorType& stdG )
	{
		int size = G.size();

		stdG = VectorType( size, 0 );
		for( int i = 0; i < size; i++ )
			stdG[i] = G( i );
	}

	/**
	 * Load vector from txt file.
	 */
	template< class ValueType >
	void Ann< ValueType >::LoadVector( VectorType& vector, const std::string& inputFileName )
	{
		std::ifstream in( inputFileName.c_str() );

		if ( !in )
		{
			std::cerr << "*** ERROR ***: Unable to open file \"" << inputFileName << "\" for input.\n";
			exit( EXIT_FAILURE );
		}

		// allocate ...
		vector.clear();

		std::string line;
		while ( getline( in, line ) )
		{
			TokType tok( line, boost::char_separator< char >( " " ) );

			// only first value in this case ...
			TokType::iterator id = tok.begin();
			vector.push_back( boost::lexical_cast< ValueType >( *id ) );
		}

		in.close();
	}

	/**
	 * Load matrix from file.
	 */
	template< class ValueType >
	void Ann< ValueType >::LoadMatrix( MatrixType& result, const std::string& matrixFileName )
	{
		typedef itk::Image< ValueType, 2 > ImageType;
		typedef itk::ImageFileReader< ImageType > ImageFileReaderType;
		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ImageConstIteratorType;

		typename ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
		reader->SetFileName( matrixFileName );
		reader->Update();
		typename ImageType::Pointer matrix = reader->GetOutput();
		typename ImageType::SizeType size = matrix->GetLargestPossibleRegion().GetSize();

		// initialize matrix ...
		result = MatrixType( size[0], VectorType( size[1], 0.0 ) );

		ImageConstIteratorType it( matrix, matrix->GetLargestPossibleRegion() );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			typename ImageType::IndexType index = it.GetIndex();
			result[index[0]][index[1]] = matrix->GetPixel( index );
		}
	}

	/**
	 * Calculates the modularity of a given clustering.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Modularity( const VectorType& clusters, const MatrixType& matrix )
	{
		if ( clusters.empty() || matrix.empty() )
			return 0.0;

		if( clusters.size() != matrix.size() )
		{
			std::cerr << "*** ERROR ***: Could not calculate modularity, ";
			std::cerr << "cluster vector size is not equal to matrix size!" << std::endl;
			exit( EXIT_FAILURE );
		}

		int maxClusters = ( *std::max_element( clusters.begin(), clusters.end() ) ) + 1;

		std::vector< ValueType > degreeIn( maxClusters );
		std::vector< ValueType > degreeOut( maxClusters );
		std::vector< ValueType > edges( maxClusters );

		ValueType m2 = 0;

		for( unsigned int i = 0; i < matrix.size(); i++ )
		{
			for( unsigned int j = 0; j < matrix[0].size(); j++ )
			{
				int c1 = clusters[i];
				int c2 = clusters[j];

				if( c1 == c2 )
					edges[c1] += matrix[i][j];

				degreeOut[c1] += matrix[i][j];
				degreeIn[c2] += matrix[i][j];
				m2 += matrix[i][j];
			}
		}

		ValueType result = 0.;
		for ( int i = 0; i < maxClusters; i++ )
		{
			result += edges[i] / m2;
			result -= ( degreeIn[i] / m2 ) * ( degreeOut[i] / m2 );
		}

		return result;
	}

	/**
	 * Return positive and negative border distances.
	 */
	template< class ValueType >
	void Ann< ValueType >::HausdorffDistance(
			ValueType& negative,
			ValueType& positive,
			const LabelImagePointerType& ref,
			const LabelImagePointerType& seg )
	{
		typename DirectedHausdorffDistanceImageFilterType::Pointer filter = DirectedHausdorffDistanceImageFilterType::New();
		filter->SetInput1( ref );
		filter->SetInput2( seg );
		filter->Update();
		negative = filter->GetAverageHausdorffDistance();

		filter->SetInput1( seg );
		filter->SetInput2( ref );
		filter->Update();
		positive = filter->GetAverageHausdorffDistance();
	}

	/**
	 * Return positive and negative border distances.
	 */
	template< class ValueType >
	void Ann< ValueType >::BorderDistance(
			ValueType& positive,
			ValueType& negative,
			const LabelImagePointerType& ref,
			const LabelImagePointerType& seg,
			const std::string& distanceMapFileName,
			const std::string& distanceMapInverseFileName
	)
	{
		// get distance map ...
		ImagePointerType distanceMap;

		DistanceMap( distanceMap, ref, distanceMapFileName );

		// get distance map invert ...
		ImagePointerType distanceMapInvert;
		LabelImagePointerType refInvert;
		Invert( ref, refInvert );

		DistanceMap( distanceMapInvert, refInvert, distanceMapInverseFileName );

		// get AND image ...
		LabelImagePointerType andImage;
		LabelAnd( andImage, seg, ref );

		// get falsePositiveImage and falseNegativeImage ...
		LabelImagePointerType falsePositiveImage;
		LabelImagePointerType falseNegativeImage;


		Subtract( falsePositiveImage, seg, andImage, "" );
		Subtract( falseNegativeImage, ref, andImage, "" );

		// get false positive and negative average distance ...
		positive = AverageValue( falsePositiveImage, distanceMap );
		negative = AverageValue( falseNegativeImage, distanceMapInvert );
	}

	/**
	 * Invert label.
	 */
	template< class ValueType >
	void Ann< ValueType >::Invert( const LabelImagePointerType& input, LabelImagePointerType& output )
	{
		LabelImageIteratorType it( input, input->GetLargestPossibleRegion() );

		output = LabelImageType::New();
		output->SetRegions( input->GetLargestPossibleRegion() );
		output->SetOrigin( input->GetOrigin() );
		output->SetSpacing( input->GetSpacing() );
		output->Allocate();
		output->FillBuffer( 0 );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			unsigned int inputValue = input->GetPixel( it.GetIndex() );

			if( inputValue == 0 )
			{
				output->SetPixel( it.GetIndex(), 1 );
			}
		}
	}

	/**
	 * Subtract image1 from image2.
	 */
	template< class ValueType >
	void Ann< ValueType >::Subtract( LabelImagePointerType& result,
			const LabelImagePointerType& image1,
			const LabelImagePointerType& image2,
			const std::string& outputFileName )
	{
		SubtractImageFilterType::Pointer subtractFilter =
			SubtractImageFilterType::New();

		subtractFilter->SetInput1( image1 );
		subtractFilter->SetInput2( image2 );
		subtractFilter->Update();
		result = subtractFilter->GetOutput();
		subtractFilter = 0;

		// write ...
		if ( !outputFileName.empty() )
		{
			WriteData( result, outputFileName );
		}
	}

	/**
	 * Distance map.
	 */
	template< class ValueType >
	void Ann< ValueType >::DistanceMap(
			ImagePointerType& distanceMap,
			const LabelImagePointerType& ref,
			const std::string& outputFileName )
	{
		// cast ref label images to ValueType image ...
		typename CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
		castFilter->SetInput( ref );
		castFilter->Update();

		typename SignedMaurerDistanceMapImageFilterType::Pointer distanceFilter =
			SignedMaurerDistanceMapImageFilterType::New();
		distanceFilter->SetInput( castFilter->GetOutput() );
		distanceFilter->SetBackgroundValue( 0 );
		distanceFilter->SetSquaredDistance( false );
		distanceFilter->UseImageSpacingOn();
		distanceFilter->Update();
		distanceMap = distanceFilter->GetOutput();

		castFilter = 0;
		distanceFilter = 0;

		// write ...
		if ( !outputFileName.empty() )
		{
			WriteData( distanceMap, outputFileName );
		}
	}

	/**
	 * Mean value in first label.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::AverageValue( const LabelImagePointerType& labels,
							const ImagePointerType& values )
	{
		typename LabelStatisticsImageFilterType::Pointer meanFilter =
			LabelStatisticsImageFilterType::New();

		// mean value ...
		meanFilter->SetLabelInput( labels );
		meanFilter->SetInput( values );
		meanFilter->Update();

		return meanFilter->GetMean( 1 );
	}

	/**
	 * Range of border distances.
	 */
	template< class ValueType >
	void Ann< ValueType >::BorderDistanceRange( VectorType& negative, VectorType& positive,
			const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		ImagePointerType seg;
		GetImage( seg, segFileName );

		BorderDistanceRange( negative, positive, ref, seg, minVoxels );
	}

	/**
	 * Range of border distances.
	 */
	template< class ValueType >
	void Ann< ValueType >::HausdorffDistanceRange( VectorType& average, VectorType& directed,
			const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		ImagePointerType seg;
		GetImage( seg, segFileName );

		HausdorffDistanceRange( average, directed, ref, seg, minVoxels );
	}

	/**
	 * Range of Hausdorff distance range.
	 */
	template< class ValueType >
	void Ann< ValueType >::HausdorffDistanceRange( VectorType& average, VectorType& directed, const LabelImagePointerType& ref,
			const ImagePointerType& seg, unsigned int minVoxels )
	{
		for ( ValueType i = 0; i < 0.99; i += 0.025 )
		{
			// threshold ...
			LabelImagePointerType thresImage;
			ThresholdToBinaryImage( seg, i, thresImage );

			ValueType posValue;
			ValueType negValue;

			// remove small blobs ...
			if ( minVoxels != 0 )
			{
				LabelImagePointerType minSizeImage;
				RemoveIsolatedComponents( thresImage, minSizeImage, minVoxels );

				HausdorffDistance( posValue, negValue, ref, minSizeImage );

				average.push_back( posValue );
				directed.push_back( negValue );
			}
			else
			{
				LabelImagePointerType minSizeImage;
				RemoveIsolatedComponents( thresImage, minSizeImage, 0 ); // all voxels ...

				HausdorffDistance( posValue, negValue, ref, minSizeImage );

				average.push_back( posValue );
				directed.push_back( negValue );
			}
		}
	}

	/**
	 * Border distance (positive and negative).
	 */
	template< class ValueType >
	void Ann< ValueType >::HausdorffDistance( ValueType& average, ValueType& directed,
			const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels )

	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		LabelImagePointerType seg;
		GetLabelImage( seg, segFileName );

		// remove small blobs ...
		if ( minVoxels != 0 )
		{
			LabelImagePointerType minSizeImage;
			RemoveIsolatedComponents( seg, minSizeImage, minVoxels );
			HausdorffDistance( average, directed, ref, minSizeImage );
		}
		else
		{
			HausdorffDistance( average, directed, ref, seg );
		}
	}

	/**
	 * Range of border distances.
	 */
	template< class ValueType >
	void Ann< ValueType >::BorderDistanceRange( VectorType& negative, VectorType& positive, const LabelImagePointerType& ref,
			const ImagePointerType& seg, unsigned int minVoxels )
	{
		for ( ValueType i = 0; i < 0.99; i += 0.025 )
		{
			// threshold ...
			LabelImagePointerType thresImage;
			ThresholdToBinaryImage( seg, i, thresImage );

			ValueType posValue;
			ValueType negValue;

			// remove small blobs ...
			if ( minVoxels != 0 )
			{
				LabelImagePointerType minSizeImage;
				RemoveIsolatedComponents( thresImage, minSizeImage, minVoxels );

				BorderDistance( posValue, negValue, ref, minSizeImage, "", "" );

				positive.push_back( posValue );
				negative.push_back( negValue );
			}
			else
			{
				LabelImagePointerType minSizeImage;
				RemoveIsolatedComponents( thresImage, minSizeImage, 0 ); // all voxels ...

				BorderDistance( posValue, negValue, ref, minSizeImage, "", "" );

				positive.push_back( posValue );
				negative.push_back( negValue );
			}
		}
	}

	/**
	 * Border distance (positive and negative).
	 */
	template< class ValueType >
	void Ann< ValueType >::BorderDistance( ValueType& positive, ValueType& negative,
			const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels,
			const std::string& distanceMapFileName,
			const std::string& distanceMapInverseFileName )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		LabelImagePointerType seg;
		GetLabelImage( seg, segFileName );

		// remove small blobs ...
		if ( minVoxels != 0 )
		{
			LabelImagePointerType minSizeImage;
			RemoveIsolatedComponents( seg, minSizeImage, minVoxels );
			BorderDistance( positive, negative, ref, minSizeImage, distanceMapFileName, distanceMapInverseFileName );
		}
		else
		{
			BorderDistance( positive, negative, ref, seg, distanceMapFileName, distanceMapInverseFileName );
		}
	}

	/**
	 * Add angles to data matrix.
	 */
	template< class ValueType >
	void Ann< ValueType >::FillAnglesData( MatrixType& features,
					const std::string& labelFileName, const std::string& maskFileName )
	{
		if ( labelFileName.empty() || maskFileName.empty() )
		{
			std::cerr << "Label file or mask file not specified!" << std::endl;
			std::cerr << "[Label file only needed when spatial/angle features are enabled]" << std::endl;
			exit( EXIT_FAILURE );
		}

		ImagePointerType label;
		GetImage( label, labelFileName );
		PointType cog;
		GetCog( cog, maskFileName );

		// initialize angle features ...
		MatrixType angles;

		// unique dim-permutations (e.g. for 3 Dims: 0-1, 0-2, 1-2)
		unsigned int featureNr = 0;
		for ( unsigned int i = 0; i < Dimension; i++ )
		{
			for( unsigned int j = i + 1; j < Dimension; j++ )
			{
				VectorType angle( features[0].size(), 0 );
				angles.push_back( angle );
				featureNr++;
			}
		}

		ImageIteratorType lit( label, label->GetLargestPossibleRegion() );
		typename ImageType::SpacingType spacing = label->GetSpacing();
		unsigned int voxel = 0;

		for ( lit.GoToBegin(); !lit.IsAtEnd(); ++lit )
		{
			typename ImageType::IndexType index = lit.GetIndex();

			ValueType labelPixel = label->GetPixel( index );
			if ( labelPixel != 0 )
			{
				// unique dim-permutations (e.g. for 3 Dims: 0-1, 0-2, 1-2)
				unsigned int round = 0;
				for ( unsigned int i = 0; i < Dimension; i++ )
				{
					for( unsigned int j = i + 1; j < Dimension; j++ )
					{
						angles[round][voxel] = Angle( index[i] * spacing[i] - cog[i], index[j] * spacing[j] - cog[j] );
						round++;
					}
				}
				voxel++;
			}
		}

		// add angle features to input features ...
		for ( unsigned int i = 0; i < featureNr; i++ )
		{
			features.push_back( angles[i] );
		}
	}

	/**
	 * Return angle for given x,y in degrees.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Angle( ValueType x, ValueType y )
	{
		return std::atan2( y, x ) * 180 / PI;
	}

	/**
	 * Range SI.
	 */
	template< class ValueType >
	void Ann< ValueType >::RangeSI( VectorType& si, const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		ImagePointerType seg;
		GetImage( seg, segFileName );

		for ( ValueType i = 0; i < 0.99; i += 0.01 )
		{
			// threshold ...
			LabelImagePointerType thresImage;
			ThresholdToBinaryImage( seg, i, thresImage );

			// remove small blobs ...
			if ( minVoxels != 0 )
			{
				LabelImagePointerType minSizeImage;
				RemoveIsolatedComponents( thresImage, minSizeImage, minVoxels );
				si.push_back( SI( ref, minSizeImage ) );
			} else
			{
				si.push_back( SI( ref, thresImage ) );
			}
		}
	}

	/**
	 * ROC Curve (two vectors as output -> matrix).
	 *
	 * specificity = ~ ( 1 - specificity ) !
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::ROC( VectorType& sensitivity, VectorType& specificity, VectorType& roc,
			const std::string& refFileName,
			const std::string& segFileName,
			const std::string& maskFileName,
			unsigned int minVoxels )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		ImagePointerType seg;
		GetImage( seg, segFileName );

		LabelImagePointerType mask;
		GetLabelImage( mask, maskFileName );

		for ( ValueType i = 0; i < 0.99; i += 0.01 )
		{
			// threshold ...
			LabelImagePointerType thres;
			ThresholdToBinaryImage( seg, i, thres );

			// remove small blobs ...
			if ( minVoxels != 0 )
			{
				LabelImagePointerType minSizeImage;
				RemoveIsolatedComponents( thres, minSizeImage, minVoxels );

				sensitivity.push_back( Sensitivity( ref, minSizeImage ) );
				specificity.push_back( 1.0 - Specificity( ref, minSizeImage, mask ) );
			}
			else
			{
				sensitivity.push_back( Sensitivity( ref, thres ) );
				specificity.push_back( 1.0 - Specificity( ref, thres, mask ) );
			}
		}

		// cumulative ROC ... (return AUC)
		return AUC( sensitivity, specificity, roc );
	}

	/**
	 * Fill ROC-vector with cumulative AUC, and return final AUC.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::AUC( const VectorType& sensitivity, const VectorType& specificity, VectorType& roc )
	{
		// =C0 + (B0-B1) * (A0+A1) / 2 -> gnumeric ...
		roc = VectorType( sensitivity.size(), 0 );

		for ( unsigned int i = 0; i < sensitivity.size() - 1; i++ )
		{
			ValueType deltaB = specificity[i] - specificity[i + 1];
			ValueType deltaA = sensitivity[i] + sensitivity[i + 1];

			roc[i + 1] = roc[i] + deltaB * deltaA / 2.0;
		}

		return roc[roc.size() - 1];
	}

	/**
	 * Threshold image (lower thr.).
	 */
	template< class ValueType >
	void Ann< ValueType >::ThresholdToBinaryImage( const ImagePointerType& input, ValueType threshold, LabelImagePointerType& output )
	{
		typename BinaryThresholdImageFilterType::Pointer filter = BinaryThresholdImageFilterType::New();
		filter->SetInput( input );
		filter->SetLowerThreshold( threshold );
		filter->SetUpperThreshold( 255 );
		filter->Update();
		output = filter->GetOutput();
		filter = 0;
	}

	/**
	 * Threshold image (lower thr.).
	 */
	template< class ValueType >
	void Ann< ValueType >::ThresholdToImage( const ImagePointerType& input, ValueType threshold, ImagePointerType& output )
	{
		typename ThresholdImageFilterType::Pointer filter = ThresholdImageFilterType::New();
		filter->SetInput( input );
		filter->ThresholdBelow( threshold );
		filter->SetOutsideValue( 0 );
		//filter->SetUpperThreshold( std::numeric_limits< ValueType >::max() );
		filter->Update();
		output = filter->GetOutput();
		filter = 0;
	}

	/**
	 * Crear label image from minimal components.
	 */
	template< class ValueType >
	void Ann< ValueType >::RemoveIsolatedComponents( const std::string& inputFileName, const std::string& outputFileName,
			unsigned int minVoxels )
	{
		LabelImagePointerType input;
		GetLabelImage( input, inputFileName );

		LabelImagePointerType output;
		RemoveIsolatedComponents( input, output, minVoxels );

		WriteData( output, outputFileName );
	}

	/**
	 * Remove isolated components.
	 */
	template< class ValueType >
	void Ann< ValueType >::RemoveIsolatedComponents( const LabelImagePointerType& input, LabelImagePointerType& output,
			unsigned int minVoxels )
	{
		// connected filter ...
		typename ConnectedComponentImageFilterType::Pointer componentFilter = ConnectedComponentImageFilterType::New();
		componentFilter->SetInput( input );
		componentFilter->Update();

		// relabel filter ...
		typename RelabelComponentImageFilterType::Pointer relabel = RelabelComponentImageFilterType::New();
		relabel->SetInput( componentFilter->GetOutput() );
		relabel->SetMinimumObjectSize( minVoxels );
		relabel->Update();

		// output ...
		output = relabel->GetOutput();

		componentFilter = 0;
		relabel = 0;
	}

	/**
	 * Write label image to file.
	 */
	template< class ValueType >
	void Ann< ValueType >::WriteData( const LabelImagePointerType& output, const std::string& outputFileName )
	{
		typename LabelImageWriterType::Pointer writer = LabelImageWriterType::New();
		writer->SetInput( output );
		writer->SetFileName( outputFileName );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not write image to: " << outputFileName << "!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Write image to file.
	 */
	template< class ValueType >
	void Ann< ValueType >::WriteData( const ImagePointerType& output, const std::string& outputFileName )
	{
		typename ImageWriterType::Pointer writer = ImageWriterType::New();
		writer->SetInput( output );
		writer->SetFileName( outputFileName );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not write image to: " << outputFileName << "!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Sensitivity.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Sensitivity( const std::string& refFileName, const std::string& segFileName )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		LabelImagePointerType seg;
		GetLabelImage( seg, segFileName );

		return Sensitivity( ref, seg );
	}

	/**
	 * Sensitivity.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Sensitivity( const LabelImagePointerType& ref, const LabelImagePointerType& seg )
	{
		ValueType tp = TP( ref, seg );
		ValueType fn = FN( ref, seg );

		if ( ( tp + fn ) != 0 )
		{
			return tp / ( tp + fn );
		} else
		{
			return 0;
		}
	}

	/**
	 * Specificity.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Specificity( const std::string& refFileName, const std::string& segFileName,
			const std::string& maskFileName )
	{
		typename LabelImageType::Pointer ref;
		GetLabelImage( ref, refFileName );

		typename LabelImageType::Pointer seg;
		GetLabelImage( seg, segFileName );

		typename LabelImageType::Pointer mask;
		GetLabelImage( mask, maskFileName );

		return Specificity( ref, seg, mask );
	}

	/**
	 * Specificity.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Specificity( const LabelImagePointerType& ref, const LabelImagePointerType& seg,
			const LabelImagePointerType& mask )
	{
		ValueType tn = TN( ref, seg, mask );
		ValueType fp = FP( ref, seg );

		if ( ( fp + tn ) != 0 )
		{
			return tn / ( fp + tn );
		} else
		{
			return 0;
		}
	}

	/**
	 * SI.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::SI( const std::string& refFileName, const std::string& segFileName )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		LabelImagePointerType seg;
		GetLabelImage( seg, segFileName );

		return SI( ref, seg );
	}

	/**
	 * SI.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::SI( const LabelImagePointerType& ref, const LabelImagePointerType& seg )
	{
		LabelImagePointerType andImage;

		// AND-operator.
		LabelAnd( andImage, ref, seg );

		ValueType andVoxels = static_cast< ValueType > ( GetNonZeroVoxels( andImage ) );
		ValueType refVoxels = static_cast< ValueType > ( GetNonZeroVoxels( ref ) );
		ValueType segVoxels = static_cast< ValueType > ( GetNonZeroVoxels( seg ) );

		return ( 2 * andVoxels ) / ( refVoxels + segVoxels );
	}

	/**
	 * PSI.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::PSI( const std::string& refFileName, const std::string& segFileName )
	{
		LabelImagePointerType ref;
		GetLabelImage( ref, refFileName );

		ImagePointerType seg;
		GetImage( seg, segFileName );

		return PSI( ref, seg );
	}

	template< class ValueType >
	ValueType Ann< ValueType >::PSI( const LabelImagePointerType& ref, const ImagePointerType& seg )
	{
		ValueType refVoxels = static_cast< ValueType > ( GetNonZeroVoxels( ref ) );
		ValueType segVoxels = static_cast< ValueType > ( GetNonZeroVoxels( seg ) );

		// AND-operator.
		ValueType andVoxels = And( ref, seg );

		return ( 2 * andVoxels ) / ( refVoxels + segVoxels );
	}

	// ########################################

	/**
	 * True positives.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::TP( const LabelImagePointerType& ref, const LabelImagePointerType& seg )
	{
		LabelImageIteratorType rit( ref, ref->GetLargestPossibleRegion() );
		LabelImageIteratorType sit( seg, seg->GetLargestPossibleRegion() );

		ValueType numerator = 0;
		ValueType denominator = 0;

		for ( rit.GoToBegin(), sit.GoToBegin(); !rit.IsAtEnd(), !sit.IsAtEnd(); ++rit, ++sit )
		{
			unsigned int refValue = ref->GetPixel( rit.GetIndex() );
			unsigned int segValue = seg->GetPixel( sit.GetIndex() );

			if ( segValue != 0 && refValue != 0 )
			{
				numerator++;
			}

			if ( refValue != 0 )
			{
				denominator++;
			}
		}

		if ( denominator == 0 )
		{
			return 0;
		} else
		{
			return numerator / denominator;
		}
	}

	/**
	 * True negatives.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::TN( const LabelImagePointerType& ref, const LabelImagePointerType& seg, const LabelImagePointerType& mask )
	{
		LabelImageIteratorType rit( ref, ref->GetLargestPossibleRegion() );
		LabelImageIteratorType sit( seg, seg->GetLargestPossibleRegion() );
		LabelImageIteratorType mit( mask, mask->GetLargestPossibleRegion() );

		ValueType numerator = 0;
		ValueType denominator = 0;

		for ( mit.GoToBegin(), rit.GoToBegin(), sit.GoToBegin(); !rit.IsAtEnd(), !sit.IsAtEnd(), !mit.IsAtEnd(); ++rit, ++sit, ++mit )
		{
			unsigned int refValue = ref->GetPixel( rit.GetIndex() );
			unsigned int segValue = seg->GetPixel( sit.GetIndex() );
			unsigned int maskValue = mask->GetPixel( mit.GetIndex() );

			if ( maskValue != 0 && refValue == 0 )
			{
				numerator++;
			}

			if ( maskValue != 0 && ( refValue == 0 || segValue == 0 ) )
			{
				denominator++;
			}

		}

		if ( denominator == 0 )
		{
			return 0;
		} else
		{
			return numerator / denominator;
		}
	}

	/**
	 * False negatives.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::FN( const LabelImagePointerType& ref, const LabelImagePointerType& seg )
	{
		LabelImageIteratorType rit( ref, ref->GetLargestPossibleRegion() );
		LabelImageIteratorType sit( seg, seg->GetLargestPossibleRegion() );

		ValueType numerator = 0;
		ValueType denominator = 0;

		for ( rit.GoToBegin(), sit.GoToBegin(); !rit.IsAtEnd(), !sit.IsAtEnd(); ++rit, ++sit )
		{
			unsigned int refValue = ref->GetPixel( rit.GetIndex() );
			unsigned int segValue = seg->GetPixel( sit.GetIndex() );

			if ( segValue == 0 && refValue != 0 )
			{
				numerator++;
			}

			if ( refValue != 0 )
			{
				denominator++;
			}
		}

		if ( denominator == 0 )
		{
			return 0;
		} else
		{
			return numerator / denominator;
		}
	}

	/**
	 * False positives.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::FP( const LabelImagePointerType& ref, const LabelImagePointerType& seg )
	{
		LabelImageIteratorType rit( ref, ref->GetLargestPossibleRegion() );
		LabelImageIteratorType sit( seg, seg->GetLargestPossibleRegion() );

		ValueType numerator = 0;
		ValueType denominator = 0;

		for ( rit.GoToBegin(), sit.GoToBegin(); !rit.IsAtEnd(), !sit.IsAtEnd(); ++rit, ++sit )
		{
			unsigned int refValue = ref->GetPixel( rit.GetIndex() );
			unsigned int segValue = seg->GetPixel( sit.GetIndex() );

			if ( segValue != 0 && refValue == 0 )
			{
				numerator++;
			}

			if ( refValue != 0 )
			{
				denominator++;
			}
		}

		if ( denominator == 0 )
		{
			return 0;
		} else
		{
			return numerator / denominator;
		}
	}

	/**
	 * Return number of nonzero voxels.
	 */
	template< class ValueType >
	unsigned int Ann< ValueType >::GetNonZeroVoxels( const LabelImagePointerType image )
	{
		LabelImageIteratorType it( image, image->GetLargestPossibleRegion() );

		unsigned int totalVoxels = 0;
		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			if ( image->GetPixel( it.GetIndex() ) != 0 )
			{
				totalVoxels++;
			}
		}

		return totalVoxels;
	}

	/**
	 * Return number + value of nonzero voxels.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::GetNonZeroVoxels( const ImagePointerType image )
	{
		ImageIteratorType it( image, image->GetLargestPossibleRegion() );

		ValueType totalVoxels = 0;
		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			ValueType value = image->GetPixel( it.GetIndex() );
			if ( value != 0 )
			{
				totalVoxels += value;
			}
		}

		return totalVoxels;
	}

	/**
	 * And operator binary images.
	 */
	template< class ValueType >
	void Ann< ValueType >::LabelAnd( LabelImageType::Pointer& result, const LabelImagePointerType image1,
			const LabelImageType::Pointer image2 )
	{
		LabelAndImageFilterType::Pointer filter = LabelAndImageFilterType::New();
		filter->SetInput1( image1 );
		filter->SetInput2( image2 );
		filter->Update();
		result = filter->GetOutput();
		filter = 0;
	}

	/**
	 * Or operator binary images.
	 */
	template< class ValueType >
	void Ann< ValueType >::LabelOr( LabelImageType::Pointer& result, const LabelImagePointerType image1,
			const LabelImageType::Pointer image2 )
	{
		LabelOrImageFilterType::Pointer filter = LabelOrImageFilterType::New();
		filter->SetInput1( image1 );
		filter->SetInput2( image2 );
		filter->Update();
		result = filter->GetOutput();
		filter = 0;
	}

	/**
	 * Return AND of images.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::And( const LabelImagePointerType ref, const ImagePointerType seg )
	{
		LabelImageIteratorType rit( ref, ref->GetLargestPossibleRegion() );
		ImageIteratorType sit( seg, seg->GetLargestPossibleRegion() );

		ValueType totalVoxels = 0;

		for ( rit.GoToBegin(), sit.GoToBegin(); !rit.IsAtEnd(), !sit.IsAtEnd(); ++rit, ++sit )
		{
			// for non-zero voxels in reference ...
			if ( ref->GetPixel( rit.GetIndex() ) != 0 )
			{
				// add probability in seg ...
				totalVoxels += seg->GetPixel( sit.GetIndex() );
			}
		}

		return totalVoxels;
	}

	/**
	 * Return image pointer.
	 */
	template< class ValueType >
	void Ann< ValueType >::GetImage( ImagePointerType& result, const std::string& imageFileName )
	{
		typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
		imageReader->SetFileName( imageFileName.c_str() );
		imageReader->Update();
		result = imageReader->GetOutput();
		imageReader = 0;
	}

	/**
	 * Return image pointer.
	 */
	template< class ValueType >
	void Ann< ValueType >::GetLabelImage( LabelImageType::Pointer& result, const std::string& imageFileName )
	{
		LabelImageReaderType::Pointer imageReader = LabelImageReaderType::New();
		imageReader->SetFileName( imageFileName.c_str() );
		imageReader->Update();
		result = imageReader->GetOutput();
		imageReader = 0;
	}

	/**
	 * Check Dimensions.
	 */
	template< class ValueType >
	void Ann< ValueType >::CheckDimensions( const std::string& image )
	{
		if ( Dimension != GetImageDimensions( image ) )
		{
			std::cerr << "Dimension of input image: " << image << " is != " << Dimension << "!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Get image dimensions.
	 */
	template< class ValueType >
	unsigned int Ann< ValueType >::GetImageDimensions( const std::string& inputFileName )
	{
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io )
		{
			std::cerr << "Could not create a valid ImageIO for: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else
		{
			io->SetFileName( inputFileName );
			io->ReadImageInformation();
			return io->GetNumberOfDimensions();
		}
	}

	/**
	 * Check Dimensions.
	 */
	template< class ValueType >
	void Ann< ValueType >::CheckDimensions( const std::vector< std::string >& images )
	{
		for ( unsigned int i = 0; i < images.size(); i++ )
		{
			CheckDimensions( images[i] );
		}
	}

	/**
	 * Read images and check for dimensions, etc.
	 */
	template< class ValueType >
	void Ann< ValueType >::ReadImages( std::vector< ImagePointerType >& images, const std::vector< std::string >& inputFileNames )
	{
		// First check if all input has the same dimension.

		CheckDimensions( inputFileNames );

		// read images ...
		for ( unsigned int i = 0; i < inputFileNames.size(); i++ )
		{
			ImagePointerType image;
			ReadImage( image, inputFileNames[i].c_str() );
			images.push_back( image );
		}
	}

	/**
	 * Read image from file.
	 */
	template< class ValueType >
	void Ann< ValueType >::ReadImage( ImagePointerType& image, const std::string& imageFileName )
	{
		CheckDimensions( imageFileName );
		typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
		imageReader->SetFileName( imageFileName.c_str() );
		imageReader->Update();
		image = imageReader->GetOutput();
		imageReader = 0;
	}

	/**
	 * Read label image from file.
	 */
	template< class ValueType >
	void Ann< ValueType >::ReadImage( LabelImagePointerType& image, const std::string& imageFileName )
	{
		CheckDimensions( imageFileName );
		typename LabelImageReaderType::Pointer imageReader = LabelImageReaderType::New();
		imageReader->SetFileName( imageFileName.c_str() );
		imageReader->Update();
		image = imageReader->GetOutput();
		imageReader = 0;
	}

	/**
	 * Return center-of-gravity for 3D volume.
	 */
	template< class ValueType >
	void Ann< ValueType >::GetCog( PointType& point, const std::string& imageFileName )
	{
		ImagePointerType image;
		GetImage( image, imageFileName );

		PointType pos;

		ValueType total_masses = 0;

		ImageIteratorType it( image, image->GetLargestPossibleRegion() );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			typename ImageType::IndexType index = it.GetIndex();
			ValueType mass = image->GetPixel( index );

			for ( unsigned int i = 0; i < Dimension; i++ )
			{
				pos[i] += mass * index[i];
			}

			total_masses += mass;
		}

		typename ImageType::SpacingType spacing = image->GetSpacing();

		for ( unsigned int i = 0; i < Dimension; i++ )
		{
			point[i] = ( pos[i] * spacing[i] ) / total_masses;
		}
	}

	/**
	 * Return matrix with spatial x,y,z.
	 */
	template< class ValueType >
	void Ann< ValueType >::FillSpatialData( MatrixType& features, const std::string& labelFileName, const std::string& maskFileName )
	{
		if ( labelFileName.empty() || maskFileName.empty() )
		{
			std::cerr << "Label file or mask file not specified!" << std::endl;
			std::cerr << "[Label file only needed when spatial features are enabled]" << std::endl;
			exit( EXIT_FAILURE );
		}

		ImagePointerType label;
		GetImage( label, labelFileName );
		PointType cog;
		GetCog( cog, maskFileName );

		// initialize distance features ...
		MatrixType spatials;

		for ( unsigned int i = 0; i < Dimension; i++ )
		{
			VectorType spatial( features[0].size(), 0 );
			spatials.push_back( spatial );
		}

		ImageIteratorType lit( label, label->GetLargestPossibleRegion() );
		typename ImageType::SpacingType spacing = label->GetSpacing();
		unsigned int voxel = 0;
		for ( lit.GoToBegin(); !lit.IsAtEnd(); ++lit )
		{
			typename ImageType::IndexType index = lit.GetIndex();

			ValueType labelPixel = label->GetPixel( index );
			if ( labelPixel != 0 )
			{
				for ( unsigned int i = 0; i < Dimension; i++ )
				{
					ValueType diff = static_cast< ValueType > ( index[i] * spacing[i] ) - static_cast< ValueType > ( cog[i] );

					spatials[i][voxel] = diff;
				}
				voxel++;
			}
		}

		// add spatial features to input features ...
		for ( unsigned int i = 0; i < Dimension; i++ )
		{
			features.push_back( spatials[i] );
		}
	}

	/**
	 * Return euclidean distance, given two points.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::EuclideanDistance( const PointType& A, const PointType& B )
	{
		ValueType result = 0;
		for ( unsigned int i = 0; i < Dimension; i++ )
		{
			result += sqrt( pow( A[i] - B[i], 2 ) );
		}
		return result;
	}

	/**
	 * Return matrix with euclidean distance added as features ...
	 */
	template< class ValueType >
	void Ann< ValueType >::FillEuclideanData( MatrixType& features, const std::string& labelFileName, const std::string& maskFileName )
	{
		if ( labelFileName.empty() || maskFileName.empty() )
		{
			std::cerr << "Label file or mask file not specified!" << std::endl;
			std::cerr << "[Label file only needed when spatial/polar features are enabled]" << std::endl;
			exit( EXIT_FAILURE );
		}

		ImagePointerType label;
		GetImage( label, labelFileName );
		PointType cog;
		GetCog( cog, maskFileName );

		// allocate space ...
		VectorType eucl( features[0].size(), 0 );

		ImageIteratorType lit( label, label->GetLargestPossibleRegion() );
		typename ImageType::SpacingType spacing = label->GetSpacing();
		unsigned int voxel = 0;
		for ( lit.GoToBegin(); !lit.IsAtEnd(); ++lit )
		{
			typename ImageType::IndexType index = lit.GetIndex();

			ValueType labelPixel = label->GetPixel( index );
			if ( labelPixel != 0 )
			{
				// euclidean distance ...
				PointType voxelPoint;
				for ( unsigned int i = 0; i < Dimension; i++ )
				{
					voxelPoint[i] = static_cast< ValueType > ( index[i] * spacing[i] );
				}
				eucl[voxel] = EuclideanDistance( voxelPoint, cog );
				voxel++;
			}
		}
		features.push_back( eucl );
	}

	/**
	 * Write data to output file.
	 */
	template< class ValueType >
	void Ann< ValueType >::WriteData( const MatrixType& results, const std::string& outputFileName )
	{
		if ( !results.empty() )
		{
			if ( !results[0].empty() )
			{
				std::ofstream out( outputFileName.c_str() );
				if ( out.fail() )
				{
					std::cerr << "*** ERROR ***: Not able to write to: " << outputFileName << "!" << std::endl;
					exit( EXIT_FAILURE );
				}

				for ( unsigned int r = 0; r < results[0].size(); r++ )
				{
					for ( unsigned c = 0; c < results.size(); c++ )
					{
						if ( c == results.size() - 1 )
						{
							out << results[c][r];
						} else
						{
							out << results[c][r] << "\t";
						}
					}
					out << std::endl;
				}

				out.close();
			}
		}
	}

	/**
	 * Write data to output file.
	 */
	template< class ValueType >
	void Ann< ValueType >::WriteData( const VectorType& results, const std::string& outputFileName )
	{
		if ( !results.empty() )
		{
			std::ofstream out( outputFileName.c_str() );
			if ( out.fail() )
			{
				std::cerr << "*** ERROR ***: Not able to write to: " << outputFileName << "!" << std::endl;
				exit( EXIT_FAILURE );
			}

			for ( unsigned i = 0; i < results.size(); i++ )
			{
				out << results[i] << std::endl;
			}

			out.close();
		}
	}

	/**
	 * Add intensity data to matrix.
	 */
	template< class ValueType >
	void Ann< ValueType >::FillIntensityData( MatrixType& results, const std::vector< ImagePointerType >& images,
			const LabelImagePointerType& labels )
	{
		LabelImageType::RegionType labelRegion = labels->GetLargestPossibleRegion();

		// Get image intensities ...
		for ( unsigned int i = 0; i < images.size(); i++ )
		{
			std::vector< ValueType > result;

			LabelImageIteratorType lit( labels, labelRegion );
			ImageIteratorType iit( images[i], images[i]->GetLargestPossibleRegion() );

			for ( iit.GoToBegin(), lit.GoToBegin(); !lit.IsAtEnd(), !iit.IsAtEnd(); ++lit, ++iit )
			{
				LabelPixelType labelPixel = labels->GetPixel( lit.GetIndex() );

				if ( labelPixel != 0 )
				{
					result.push_back( images[i]->GetPixel( lit.GetIndex() ) );
				}
			}

			results.push_back( result );
		}
	}

	/**
	 * Fill intensity data to matrix from given input files.
	 */
	template< class ValueType >
	void Ann< ValueType >::FillIntensityData( MatrixType& features, const std::vector< std::string >& inputFileNames,
			const std::string& labelFileName )
	{

		std::vector< ImagePointerType > images;
		ReadImages( images, inputFileNames );
		LabelImagePointerType label;
		ReadImage( label, labelFileName );

		FillIntensityData( features, images, label );
	}

	/**
	 * Add class labels to input matrix.
	 */
	template< class ValueType >
	void Ann< ValueType >::FillClassLabels( MatrixType& results, const LabelImagePointerType& labels )
	{
		LabelImageType::RegionType labelRegion = labels->GetLargestPossibleRegion();

		// Get label values ...
		std::vector< ValueType > labelIndices;

		LabelImageIteratorType lit( labels, labelRegion );
		for ( lit.GoToBegin(); !lit.IsAtEnd(); ++lit )
		{
			LabelPixelType labelPixel = labels->GetPixel( lit.GetIndex() );

			if ( labelPixel != 0 )
			{
				labelIndices.push_back( labelPixel );
			}
		}
		results.push_back( labelIndices );
	}

	/**
	 * Normalize all matrix vectors (columns).
	 * N = X - mean / sd.
	 */
	template< class ValueType >
	void Ann< ValueType >::Normalize( MatrixType& dataMatrix )
	{
		for ( unsigned int i = 0; i < dataMatrix.size(); i++ )
		{
			Normalize( dataMatrix[i] );
		}
	}

	/**
	 * Normalize vector.
	 * N = X - mean / sd.
	 */
	template< class ValueType >
	void Ann< ValueType >::Normalize( VectorType& V )
	{
		ValueType mean = std::accumulate( V.begin(), V.end(), static_cast< ValueType > ( 0 ) );
		mean /= V.size();
		ValueType sd = GetSD( V );

		std::transform( V.begin(), V.end(), V.begin(), std::bind2nd( std::minus< ValueType >(), mean ) );
		std::transform( V.begin(), V.end(), V.begin(), std::bind2nd( std::divides< ValueType >(), sd ) );
	}

	/**
	 * Return standard deviation.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::GetSD( VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::variance > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
		{
			acc( V[i] );
		}
		return sqrt( variance( acc ) );
	}

	/**
	 * Write matrix to images.
	 */
	template< class ValueType >
	void Ann< ValueType >::Project2Image( const MatrixType& M, const std::string& maskFileName, const std::string outputFileName,
			const std::vector< int >& ids )
	{
		// single output ...
		if ( ( ids.empty() ) && ( M.size() == 1 ) )
		{
			Project2Image( M[0], maskFileName, outputFileName );
			return;
		}

		if ( M.size() != ids.size() )
		{
			std::cerr << "*** ERROR ***: Could not write images, "
				"id vector size does not match with matrix size!" << std::endl;
		}

		for ( unsigned int i = 0; i < M.size(); i++ )
		{
			std::stringstream ss;
			// run twice in case .nii.gz is specified (double extension ... )
			std::string root = boost::filesystem::path( outputFileName ).replace_extension( "" ).string();
			root = boost::filesystem::path( root ).replace_extension( "" ).string();

			ss << root;
			ss << "_" << ids[i] << ".nii.gz";

			Project2Image( M[i], maskFileName, ss.str() );
		}
	}

	/**
	 * Project vector onto image.
	 */
	template< class ValueType >
	void Ann< ValueType >::Project2Image( const VectorType& V, const std::string& maskFileName, const std::string& outputFileName )
	{
		typename ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName( maskFileName.c_str() );
		reader->Update();
		ImagePointerType mask = reader->GetOutput();
		reader = 0;
		typename ImageType::RegionType maskRegion = mask->GetLargestPossibleRegion();
		ImagePointerType output = ImageType::New();

		output->CopyInformation( mask );
		output->SetRegions( maskRegion );
		output->Allocate();
		output->FillBuffer( 0 );

		unsigned int i = 0;
		ImageIteratorType mit( mask, maskRegion );
		for ( mit.GoToBegin(); !mit.IsAtEnd(); ++mit )
		{
			if ( mask->GetPixel( mit.GetIndex() ) != 0 )
			{
				output->SetPixel( mit.GetIndex(), V[i] );
				i++;
			}
		}

		typename ImageWriterType::Pointer writer = ImageWriterType::New();
		writer->SetFileName( outputFileName.c_str() );
		writer->SetInput( output );
		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: Could not write output to: " << outputFileName << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Read training data into matrix (data) and vector(classes).
	 */
	template< class ValueType >
	void Ann< ValueType >::ReadTrainingData( const std::string& trainingDataFileName, MatrixType& trainingData,
			VectorType& trainingClasses, double percentage )
	{
		std::string sep = "\t";

		unsigned int dims = GetNumberOfTrainingDims( trainingDataFileName );
		unsigned int points = GetNumberOfTrainingPoints( trainingDataFileName );

		// initialize matrix ...
		for ( unsigned int i = 1; i < dims; i++ )
		{
			std::vector< ValueType > V( points, 0 );
			trainingData.push_back( V );
		}

		// initialize vector ...
		trainingClasses.resize( points, 0 );

		// open training data ...
		std::ifstream in( trainingDataFileName.c_str() );
		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read labels from: " << trainingDataFileName << std::endl;
			exit( EXIT_FAILURE );
		}

		// insert all values into matrix, except the last one (classes id) which is inserted in
		// the trainingClasses vector ...
		std::string line;
		unsigned int rowIndex = 0;
		while ( getline( in, line ) )
		{
			TokType tok( line, boost::char_separator< char >( sep.c_str() ) );
			try
			{
				unsigned int colIndex = 0;
				for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
				{
					if ( colIndex == dims - 1 ) // last column ...
					{
						trainingClasses[rowIndex] = boost::lexical_cast< ValueType >( *id );
					} else
					{
						//trainingData[colIndex - 1][rowIndex] = boost::lexical_cast< ValueType >( *id );
						trainingData[colIndex][rowIndex] = boost::lexical_cast< ValueType >( *id );
					}
					colIndex++;
				}
				rowIndex++;
			} catch ( boost::bad_lexical_cast& e )
			{
				std::cout << "*** WARNING ***: bad lexical cast during training data parsing!" << std::endl;
			}
		}
		in.close();

		if ( ( percentage > 0 ) && ( percentage < 100 ) )
		{
			ResizeTrainingData( trainingData, trainingClasses, percentage );
		}
	}

	/**
	 * Reduce number of training-points to given percentage.
	 */
	template< class ValueType >
	void Ann< ValueType >::ResizeTrainingData( MatrixType& trainingData, VectorType& trainingClasses, double percentage )
	{
		unsigned int pointsToKeep = ( trainingClasses.size() * percentage ) / 100.0;

		// initialize vector ...
		VectorType subClasses = VectorType( pointsToKeep, 0 );

		// initialize matrix ...
		MatrixType subData;
		for ( unsigned int i = 0; i < trainingData.size(); i++ )
		{
			VectorType tmp( pointsToKeep, 0 );
			subData.push_back( tmp );
		}

		// create generator ...
		random_number_type generator( time( 0 ) );

		VectorType indices;
		GetRandomValues( trainingClasses, generator, indices, pointsToKeep );

		// insert random samples ...
		for ( unsigned int i = 0; i < indices.size(); i++ )
		{
			subClasses[i] = trainingClasses[indices[i]];

			for ( unsigned int j = 0; j < subData.size(); j++ )
			{
				subData[j][i] = trainingData[j][indices[i]];
			}
		}

		trainingData = subData;
		trainingClasses = subClasses;
	}

	/**
	 * Return uniform random values vector from inputs, with similar size.
	 */
	template< class ValueType >
	void Ann< ValueType >::GetRandomValues( const VectorType& inputs, random_number_type& generator, VectorType& results,
			unsigned int total )
	{
		// uniform distribution: 0 -> largest index ...
		int_distribution_type int_uni_dist( 0, inputs.size() - 1 );

		// generator ...
		int_generator_type int_distribution( generator, int_uni_dist );

		for ( unsigned int i = 0; i < total; i++ )
		{
			results.push_back( int_distribution() );
		}
	}

	/**
	 * Return number of training data dimensions ...
	 */
	template< class ValueType >
	unsigned int Ann< ValueType >::GetNumberOfTrainingDims( const std::string& trainingData )
	{
		std::string sep = "\t";
		std::ifstream in( trainingData.c_str() );

		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read labels from: " << trainingData << std::endl;
			exit( EXIT_FAILURE );
		}

		std::string line;
		getline( in, line );

		TokType tok( line, boost::char_separator< char >( sep.c_str() ) );

		unsigned int results = 0;
		try
		{
			for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
			{
				results++;
			}

		} catch ( boost::bad_lexical_cast& e )
		{
			std::cout << "*** WARNING ***: bad lexical cast during training data parsing!" << std::endl;
		}
		in.close();
		return results;
	}

	/**
	 * Return number of training data points ...
	 */
	template< class ValueType >
	unsigned int Ann< ValueType >::GetNumberOfTrainingPoints( const std::string& trainingData )
	{
		std::ifstream in( trainingData.c_str() );

		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read labels from: " << trainingData << std::endl;
			exit( EXIT_FAILURE );
		}
		std::string line;
		unsigned int result = 0;

		while ( getline( in, line ) )
		{
			result++;
		}

		in.close();
		return result;
	}

	/**
	 * Query all rows in queryData and return results as queryResultMatrix.
	 * If probabilistic vector is empty, only one column is filled with k-nn.
	 * If probabilistic classes are specfied, each class represents one
	 * output column in the queryResultMatrix.
	 */
	template< class ValueType >
	void Ann< ValueType >::Query( MatrixType& trainingData, VectorType& trainingClasses, MatrixType& queryData,
			MatrixType& queryResultMatrix, unsigned int k, double errorBound, std::vector< int > probabilistic,
			const std::string& dumpFileName, bool prioritySearch )
	{
		if ( trainingData.size() != queryData.size() )
		{
			std::cerr << "Feature space dimension for query data does not match with training data!" << std::endl;
			exit( EXIT_FAILURE );
		}

		ANNkd_tree* kdTree = GetKDTree( trainingData );

		unsigned int nPts = queryData[0].size();
		unsigned int dim = queryData.size();

		// initialize result matrix ...
		if ( !probabilistic.empty() )
		{
			for ( unsigned int i = 0; i < probabilistic.size(); i++ )
			{
				queryResultMatrix.push_back( VectorType( nPts, 0 ) );
			}
		} else
		{
			queryResultMatrix.push_back( VectorType( nPts, 0 ) ); // fill one column with 0 ...
		}

		// query ...
		for ( unsigned int i = 0; i < nPts; i++ )
		{
			ANNpoint queryPt = annAllocPt( dim );
			ANNidxArray nnIdx = new ANNidx[k];
			ANNdistArray dists = new ANNdist[k];

			for ( unsigned int j = 0; j < dim; j++ )
			{
				queryPt[j] = static_cast< ANNdist > ( queryData[j][i] );
			}

			kdTree->annkSearch( queryPt, k, nnIdx, dists, errorBound );

			if ( prioritySearch )
			{
				kdTree->annkPriSearch( queryPt, k, nnIdx, dists, errorBound );
			}

			if ( !probabilistic.empty() )
			{
				for ( unsigned p = 0; p < probabilistic.size(); p++ )
				{
					queryResultMatrix[p][i] = Probabilistic( trainingClasses, probabilistic[p], k, nnIdx );
				}
			} else
			{
				queryResultMatrix[0][i] = trainingClasses[nnIdx[k - 1]];
			}

			delete[] nnIdx;
			delete[] dists;
		}

		if ( !dumpFileName.empty() )
		{
			DumpTree( dumpFileName, kdTree );
		}

		// clean things up
		delete kdTree;
		annClose();
	}

	/**
	 * Dump kd-tree to file.
	 */
	template< class ValueType >
	void Ann< ValueType >::DumpTree( const std::string& dumpFileName, ANNkd_tree* kdTree )
	{
		std::ofstream out( ( dumpFileName + ".dmp" ).c_str() );
		if ( out.fail() )
		{
			std::cerr << "*** ERROR ***: Not able to write to: " << dumpFileName << "!" << std::endl;
			exit( EXIT_FAILURE );
		}

		ANNbool dumpPoints = ANNtrue;
		kdTree->Dump( dumpPoints, out );
		std::cout << "To convert the dump to " << dumpFileName << ".fig: 'ann2fig " << dumpFileName << "'" << std::endl;
	}

	/**
	 * Return probability of given voxel of being part of the classNumber region,
	 * by taking the K learning samples into account.
	 */
	template< class ValueType >
	ValueType Ann< ValueType >::Probabilistic( const VectorType& trainingClasses, unsigned int classNumber, unsigned int k,
			const ANNidxArray& nnIdx )
	{
		/**
		 * For each k check if classNumber is specified ...
		 */
		unsigned int total = 0;
		for ( unsigned int i = 0; i < k - 1; i++ )
		{
			if ( classNumber == trainingClasses[nnIdx[i]] )
			{
				total++;
			}
		}

		if ( total != 0 )
		{
			return ( static_cast< ValueType > ( total ) / static_cast< ValueType > ( k - 1 ) );
		} else
		{
			return 0;
		}
	}

	/**
	 * Convert training data into a KD-tree.
	 */
	template< class ValueType >
	ANNkd_tree* Ann< ValueType >::GetKDTree( const MatrixType& trainingData )
	{
		unsigned int nPts = trainingData[0].size();
		unsigned int dim = trainingData.size();

		ANNpointArray dataPts = annAllocPts( nPts, dim ); // allocate data points

		for ( unsigned int i = 0; i < nPts; i++ )
		{
			for ( unsigned int j = 0; j < dim; j++ )
			{
				dataPts[i][j] = static_cast< double > ( trainingData[j][i] );
			}
		}

		return new ANNkd_tree( dataPts, nPts, dim );
	}
} // end namespace ann

template class ann::Ann< float >;
template class ann::Ann< double >;
