#include "graphCommon.h"

namespace graph
{
	/**
	 * Return max value.
	 */
	template< class ValueType >
	unsigned int Graph< ValueType >::GetMaxValue4DImage( const std::string& inputFileName )
	{
		typename MinMaxCalcType::Pointer filter = MinMaxCalcType::New();
		Image4DPointerType input;
		Load4DImage( input, inputFileName );
		filter->SetInput( input );
		filter->Update();
		return filter->GetMaximum();
	}

	/**
	 * Load 3d image.
	 */
	template< class ValueType >
	void Graph< ValueType >::Load3DImage( Image3DPointerType& result, const std::string& filename )
	{
		typedef itk::ImageFileReader< Image3DType > ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();
		result = reader->GetOutput();
	}

	/**
	 * Load 4d image.
	 */
	template< class ValueType >
	void Graph< ValueType >::Load4DImage( Image4DPointerType& result, const std::string& filename )
	{
		typedef itk::ImageFileReader< Image4DType > ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename.c_str() );
		reader->Update();
		result = reader->GetOutput();
	}

	/**
	 * Write 3d image.
	 */
	template< class ValueType >
	void Graph< ValueType >::Write3DImage( const Image3DPointerType& output, const std::string& filename )
	{
		typedef itk::ImageFileWriter< Image3DType > WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( filename.c_str() );
		writer->SetInput( output );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not write " << output << std::endl;
			e.GetDescription();
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Write 4d image.
	 */
	template< class ValueType >
	void Graph< ValueType >::Write4DImage( const Image4DPointerType& output, const std::string& filename )
	{
		typedef itk::ImageFileWriter< Image4DType > WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( filename.c_str() );
		writer->SetInput( output );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not write " << output << std::endl;
			e.GetDescription();
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Return matrix with vectors containing time-series for all non-zero mask voxels.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetTimeSeries( MatrixType& matrix, const std::string& inputFileName, const std::string& maskFileName )
	{
		Image4DPointerType input;
		Load4DImage( input, inputFileName );

		Image3DPointerType mask;
		Load3DImage( mask, maskFileName );

		GetTimeSeries( matrix, input, mask );
	}

	/**
	 * Return matrix with vectors containing time-series for all non-zero mask voxels.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetTimeSeries( MatrixType& matrix, const Image4DPointerType& input, const Image3DPointerType& mask )
	{
		// Get 4D-size
		typename Image4DType::SizeType size4D = input->GetLargestPossibleRegion().GetSize();
		unsigned int timePoints = size4D[3];

		// get total number of non-zero voxels (to allocate matrix...)
		unsigned int numVoxels = 0;
		Iterator3DType it( mask, mask->GetLargestPossibleRegion() );
		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
			if ( mask->GetPixel( it.GetIndex() ) != 0 )
				numVoxels++;

		// allocate ...
		matrix = MatrixType( numVoxels, timePoints );

		// Iterate over 4D image using mask...
		Iterator3DType mit( mask, mask->GetLargestPossibleRegion() );
		unsigned index = 0;
		for ( mit.GoToBegin(); !mit.IsAtEnd(); ++mit )
		{
			if ( mask->GetPixel( mit.GetIndex() ) != 0 )
			{
				typename Image3DType::IndexType index3D = mit.GetIndex();
				typename Image4DType::IndexType index4D;
				index4D[0] = index3D[0];
				index4D[1] = index3D[1];
				index4D[2] = index3D[2];

				for ( unsigned int i = 0; i < timePoints; i++ )
				{
					index4D[3] = i;
					matrix[index][i] =input->GetPixel( index4D );
				}
				index++;
			}
		}
	}

	/**
	 * Insert values of a given 3D image into a 4D image at given position.
	 */
	template< class ValueType >
	void Graph< ValueType >::Insert3DImageIn4DImage( Image4DPointerType& image4D, const Image3DPointerType& image3D, unsigned int position )
	{
		// Get 4D-size
		typename Image4DType::SizeType size4D = image4D->GetLargestPossibleRegion().GetSize();
		unsigned int volumes = size4D[3];

		// size check ...
		if( position > volumes )
		{
			std::cerr << "*** Not able to insert 3D image into 4D image: out of bounds!" << std::endl;
			exit( EXIT_FAILURE );
		}

		// Iterate over 3D image ...
		Iterator3DType it( image3D, image3D->GetLargestPossibleRegion() );

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			typename Image3DType::IndexType index3D = it.GetIndex();
			typename Image4DType::IndexType index4D;
			index4D[0] = index3D[0];
			index4D[1] = index3D[1];
			index4D[2] = index3D[2];
			index4D[3] = position;

			image4D->SetPixel( index4D, image3D->GetPixel( index3D ) );
		}
	}

	/**
	 * Return empty image from given input file image.
	 */
	template< class ValueType >
	void Graph< ValueType >::Zero3DImage( Image3DPointerType& output, const std::string& referenceFileName )
	{
		Image3DPointerType reference;
		Load3DImage( reference, referenceFileName );
		output = Image3DType::New();
		output->CopyInformation( reference );
		output->SetBufferedRegion( reference->GetBufferedRegion() );
		output->SetRequestedRegion( reference->GetRequestedRegion() );
		output->Allocate();
		output->FillBuffer( 0 );
	}

	/**
	 * Return empty 4D image from given 3D input file image and size (z).
	 */
	template< class ValueType >
	void Graph< ValueType >::Zero4DImage( Image4DPointerType& output, const std::string& referenceFileName,
			unsigned int size4thDimension )
	{
		Image3DPointerType reference;
		Load3DImage( reference, referenceFileName );
		typename Image3DType::SizeType size3D = reference->GetLargestPossibleRegion().GetSize();
		typename Image3DType::SpacingType spacing3D = reference->GetSpacing();
		typename Image3DType::PointType origin3D = reference->GetOrigin();

		output = Image4DType::New();
		typename Image4DType::RegionType region4D;
		typename Image4DType::SizeType size4D;
		typename Image4DType::IndexType index4D;
		size4D[0] = size3D[0];
		size4D[1] = size3D[1];
		size4D[2] = size3D[2];
		size4D[3] = size4thDimension;

		index4D[0] = 0;
		index4D[1] = 0;
		index4D[2] = 0;
		index4D[3] = 0;

		region4D.SetSize( size4D );
		region4D.SetIndex( index4D );

		typename Image4DType::SpacingType spacing4D;
		spacing4D[0] = spacing3D[0];
		spacing4D[1] = spacing3D[1];
		spacing4D[2] = spacing3D[2];
		spacing4D[3] = 1;

		typename Image4DType::PointType origin4D;
		origin4D[0] = origin3D[0];
		origin4D[1] = origin3D[1];
		origin4D[2] = origin3D[2];

		output->SetRegions( region4D );
		output->SetSpacing( spacing4D );
		output->SetOrigin( origin4D );

		output->Allocate();
		output->FillBuffer( 0 );
	}

	/**
	 * Return image filled with values from vector and information from reference.
	 */
	template< class ValueType >
	void Graph< ValueType >::FillImage( Image3DPointerType& fill, const std::string& referenceFileName, const VectorType& values )
	{
		Image3DPointerType image;
		Load3DImage( image, referenceFileName );
		Zero3DImage( fill, referenceFileName );

		// Iterate over images...
		Iterator3DType iit( image, image->GetLargestPossibleRegion() );
		Iterator3DType fit( fill, fill->GetLargestPossibleRegion() );

		unsigned int index = 0;
		for ( iit.GoToBegin(), fit.GoToBegin(); !iit.IsAtEnd(), !fit.IsAtEnd(); ++iit, ++fit )
		{
			if ( ( image->GetPixel( iit.GetIndex() ) != 0 ) & ( index < values.size() ) )
			{
				fill->SetPixel( fit.GetIndex(), values[index] );
				++index;
			}
		}
	}

	/**
	 * Write majority vote label image to output file.
	 */
	template< class ValueType >
	void Graph< ValueType >::MajorityVote( const std::string& inputFileName,
			const std::string& maskFileName, const std::string& outputFileName )
	{
		MatrixType timeSeries;
		GetTimeSeries( timeSeries, inputFileName, maskFileName );

		unsigned int nLabels = GetMaxValue4DImage( inputFileName );
		VectorType vote( timeSeries.rows() );
		for ( unsigned int i = 0; i < timeSeries.rows(); i++ )
		{
			VectorType V = timeSeries.get_row( i );
			vote[i] = MajorityVote( V, nLabels );
		}

		Image3DPointerType fill;
		FillImage( fill, maskFileName, vote );
		Write3DImage( fill, outputFileName );
	}

	/**
	 * Return most prominent label (if equal first is returned).
	 */
	template< class ValueType >
	unsigned int Graph< ValueType >::MajorityVote( const VectorType& a, unsigned int nLabels )
	{
		if ( a.empty() )
		{
			std::cerr << "*** ERROR ***: could not calculate majority vote (vector empty)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		unsigned int label = 0;
		unsigned int minSize = 0;

		for( unsigned int i = 1; i <= nLabels; i++ )
		{
			unsigned int min = 0;
			for( unsigned int j = 0; j < a.size(); j++ )
			{
				if( a[j] == i )
					min++;
			}

			if( min > minSize )
			{
				minSize = min;
				label = i;
			}
		}

		return label;
	}

	/**
	 * Discrete entropy.
	 */
	template< class ValueType >
	void Graph< ValueType >::Entropy( const std::string& inputFileName,
			const std::string& maskFileName, const std::string& outputFileName )
	{
		MatrixType timeSeries;
		GetTimeSeries( timeSeries, inputFileName, maskFileName );
		unsigned int nLabels = GetMaxValue4DImage( inputFileName );
		VectorType entropy( timeSeries.rows() );
		Entropy( entropy, timeSeries, nLabels );

		Image3DPointerType fill;
		FillImage( fill, maskFileName, entropy );
		Write3DImage( fill, outputFileName );
	}

	/**
	 * Discrete entropy.
	 */
	template< class ValueType >
	void Graph< ValueType >::Entropy( VectorType& entropy, const MatrixType& timeSeries, unsigned int nLabels )
	{
		entropy = VectorType( timeSeries.rows() );
		for ( unsigned int i = 0; i < timeSeries.rows(); i++ )
		{
			VectorType V = timeSeries.get_row( i );
			entropy[i] = Entropy( V, nLabels );
		}
	}

	/**
	 * Return discrete entropy.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::Entropy( const VectorType& a, unsigned int nLabels )
	{
		if ( a.empty() )
		{
			std::cerr << "*** ERROR ***: could not calculate discrete entropy (vector empty)!" << std::endl;
			exit( EXIT_FAILURE );
		}

		ValueType H = 0;
		for( unsigned int i = 1; i <= nLabels; i++ )
		{
			ValueType prob = 0;
			for( unsigned int j = 0; j < a.size(); j++ )
			{
				if( a[j] == i )
					prob += 1.;
			}

			prob /= static_cast< ValueType >( a.size() );

			if( prob != 0 )
				H += prob * std::log( prob );
		}

		return -H;
	}

	/**
	 * Run Kmeans++.
	 *
	 * data = full matrix.
	 * N = number of clusters.
	 */
	template< class ValueType >
	void Graph< ValueType >::KMeansPP( const MatrixType& data, int numberOfClusters, VectorType& clusters )
	{
		int N = data.rows();
		int dim = numberOfClusters;
		gsis::Vector ** x = gsis::Kmeans::readInputData( data, N, dim );
		gsis::Kmeans kmeans( numberOfClusters, N, dim, x );
		kmeans.setThreshold( exp( -10 ) );
		kmeans.setKmeansPP( true );
		kmeans.clustering();
		clusters = kmeans.getResults();
	}

	/**
	 * Write eigenfiles. Ignore first eigenvector!
	 */
	template< class ValueType >
	void Graph< ValueType >::WriteEigenFiles( const std::string& outputRoot, unsigned int numberOfEigenVectors,
			const MatrixType& eigenVectors, const VectorType& eigenValues )
	{
		WriteVectorToFile( outputRoot + "_values.txt", eigenValues );

		if( numberOfEigenVectors > eigenVectors.size() )
		{
			std::cerr << "*** ERROR ***: Number of eigenvectors to write too large -> Resetted to maximum" << std::endl;
			numberOfEigenVectors = eigenVectors.size();
		}

		for( unsigned int i = 1; i <= numberOfEigenVectors; i++ ) // ignore first eigenvector.
		{
			std::stringstream ss;
			ss << outputRoot << "_vector_" << i << ".txt";
			WriteVectorToFile( ss.str(), eigenVectors.get_column( i ) );
		}
	}

	/**
	 * Calculate laplacian.
	 */
	template< class ValueType >
	void Graph< ValueType >::CalculateLaplacian( MatrixType& L, const MatrixType& G, DiagMatrixType& D, bool normalize )
	{
		int rows = G.rows();
		int cols = G.cols();

		L = MatrixType( rows, cols );

		// calculate node degrees
		VectorType degrees( rows );
		for ( int i = 0; i < rows; ++i )
		{
			ValueType d = 0;
			for ( int j = 0; j < cols; ++j )
			{
				d += G( i, j );
			}
			degrees[ i ] = d;
		}

		// construct diagonal degree matrix.
		D = DiagMatrixType( degrees );

		if ( !normalize )
		{
			// D: diagonal degrees matrix
			// A: adjacency matrix
			// L: Laplacian matrix, such that:
			// L = D - A
			for ( int i = 0; i < rows; ++i )
			{
				for ( int j = 0; j < cols; ++j )
				{
					if( i == j )
						L( i, j ) = degrees[ i ];
					else
						L( i, j ) = G( i, j ) * -1.;
				}
			}
		}
		else // normalized laplacian matrix: http://mathworld.wolfram.com/LaplacianMatrix.html
		{
			for ( int i = 0; i < rows; ++i )
			{
				for ( int j = 0; j < cols; ++j )
				{
					if ( i == j )
					{
						L( i, j ) = 1.;
					}
					else
					{
						if ( (degrees[i] != 0) & (degrees[j] != 0) )
						{
							L( i, j ) = -1. / std::sqrt( degrees[i] * degrees[j] );
						}
						else
						{
							L( i, j ) = 0;
						}
					}
				}
			}
		}
	}

	/**
	 * Generalized eigen values/vectors calculation.
	 *
	 * The eigenvalues are sorted into increasing order (of value, not absolute value).
	 */
	template< class ValueType >
	void Graph< ValueType >::GeneralizedEigenCalculation(
			const MatrixType& L,
			const DiagMatrixType& DegreeMatrix,
			MatrixType& eigenVectors,
			VectorType& eigenValues )
	{
		unsigned int rows = L.rows( );
		unsigned int cols = L.cols( );

		if ( rows != cols )
		{
			std::cerr << "*** ERROR ***: Matrix is not squared!"<< std::endl;
			exit( EXIT_FAILURE );
		}
		else if( ( rows != DegreeMatrix.rows() ) | ( cols != DegreeMatrix.cols() ) )
		{
			std::cerr << "*** ERROR ***: Laplacian matrix and degree matrix do not have equal size!"<< std::endl;
			exit( EXIT_FAILURE );
		}

		// Solves the generalized eigenproblem Ax=Bx ...
		vnl_generalized_eigensystem system = vnl_generalized_eigensystem( L, DegreeMatrix );

		eigenValues = system.D.diagonal();
		eigenVectors = system.V;
	}

	/**
	 * Generalized eigen values/vectors calculation based on sparse matrix.
	 *
	 * The eigenvalues are sorted into increasing order (of value, not absolute value).
	 */
	template< class ValueType >
	void Graph< ValueType >::GeneralizedEigenCalculation(
			vnl_sparse_matrix< ValueType >& L,
			MatrixType& eigenVectors,
			VectorType& eigenValues,
			int numberOfEigenVectors )
	{
		unsigned int rows = L.rows( );
		unsigned int cols = L.cols( );

		if ( rows != cols )
		{
			std::cerr << "*** ERROR ***: Matrix is not squared!"<< std::endl;
			exit( EXIT_FAILURE );
		}

		vnl_sparse_symmetric_eigensystem eigenSystem;

		std::cout << "Starting calculating pairs..." << std::endl;

		eigenSystem.CalculateNPairs( L, numberOfEigenVectors, true, 3 );

		std::cout << "Calculation done."<< std::endl;

		eigenVectors = MatrixType( rows, numberOfEigenVectors );
		eigenValues = VectorType( numberOfEigenVectors );

		// TODO not sure if set_row should be set_column or not....
		for( int i = 0; i < numberOfEigenVectors; i++ )
		{
			eigenVectors.set_row( i, eigenSystem.get_eigenvector( i ) );
			eigenValues( i ) = eigenSystem.get_eigenvalue( i );
		}
	}

	/**
	 * Write matrix to file.
	 */
	template< class ValueType >
	void Graph< ValueType >::WriteMatrixToFile( const MatrixType& matrix,
			const std::string& outputFileName )
	{
		ImagePointerType image = ImageType::New();
		typename ImageType::RegionType region;
		typename ImageType::SizeType size;
		size[0] = matrix.rows();
		size[1] = matrix.cols();
		region.SetSize( size );

		image->SetRegions( region );
		image->Allocate();
		image->FillBuffer( 0 );

		// Construct Image object
		for( unsigned int i = 0; i < matrix.rows(); i++ )
		{
			for( unsigned int j = 0; j < matrix.cols(); j++ )
			{
				typename ImageType::IndexType index;
				index[0] = i;
				index[1] = j;
				image->SetPixel( index, matrix( i, j ) );
			}
		}

		WriteImageToFile( image, outputFileName );
	}

	/**
	 * Write image to file.
	 */
	template< class ValueType >
	void Graph< ValueType >::WriteImageToFile( const ImagePointerType& image,
			const std::string& outputFileName )
	{
		typename WriterType::Pointer writer = WriterType::New( );
		writer->SetFileName( outputFileName.c_str( ) );
		writer->SetInput( image );
		try
		{
			writer->Update( );
		}
		catch( itk::ExceptionObject& e )
		{
			std::cerr << "*** ERROR ***: could not write: " << outputFileName << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Read matrix from file.
	 */
	template< class ValueType >
	void Graph< ValueType >::ReadMatrix( ImagePointerType& matrix, const std::string& inputFileName )
	{
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( inputFileName.c_str() );
		reader->Update();
		matrix = reader->GetOutput();
	}

	/**
	 * Return MatrixType from ImagePointerType.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetMatrix( MatrixType& output, const ImagePointerType& input )
	{
		ValueType* buffer = input->GetPixelContainer()->GetBufferPointer( );
		typename ImageType::RegionType region = input->GetLargestPossibleRegion( );
		typename ImageType::SizeType size = region.GetSize();
		unsigned int rows = size[ 0 ];
		unsigned int cols = size[ 1 ];
		DataMatrixType G( rows, cols, buffer );
		output = MatrixType( rows, cols );

		for( unsigned int i = 0; i < rows; i++ )
		{
			for( unsigned int j = 0; j < cols; j++ )
			{
				output( i, j ) = G( i, j );
			}
		}
	}

	/**
	 * Return (binary) matrix.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetMatrix( const ImagePointerType& image,
			ImagePointerType& output, ValueType threshold,
			ValueType K, bool fixed, bool weighted )
	{
		ValueType* buffer = image->GetPixelContainer()->GetBufferPointer();
		typename ImageType::RegionType region = image->GetLargestPossibleRegion();
		typename ImageType::SizeType size = region.GetSize();
		int rows = size[0];
		int cols = size[1];

		DataMatrixType data( rows, cols, buffer );

		// output image ...
		output = ImageType::New();
		output->CopyInformation( image );
		output->SetRegions( image->GetLargestPossibleRegion() );
		output->Allocate();
		output->FillBuffer( 0 );

		// simple thresholding ...
		if ( !fixed )
		{
			GetMatrix( output, data, threshold, weighted );
		}
		else // fixed K ...
		{
			// get real degree ...
			VectorType degree;
			Degree( degree, data, 0.0, true, true, false, "", "" );
			ValueType mean = degree.mean();

			// check ...
			if ( K > mean )
			{
				std::cerr << "Given K is larger than maximum degree!" << std::endl;
				exit( EXIT_FAILURE );
			}

			ValueType max = GetMaxValue( data );
			ValueType steps = max / 1000.0;

			for( ValueType t = 0; t < max; t += steps )
			{
				VectorType degree;
				Degree( degree, data, t, true, true, false, "", "" );

				if ( K >= degree.mean() )
				{
					GetMatrix( output, data, t, weighted );
					//std::cout << t << std::endl;
					return;
				}
			}
		}
	}

	/**
	 * Return (binary) thresholded matrix.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetMatrix( ImagePointerType& output,
			const DataMatrixType& data, ValueType threshold,
			bool weighted )
	{
		unsigned int rows = data.rows();
		unsigned int cols = data.cols();

		for ( unsigned int i = 0; i < rows; ++i )
		{
			for ( unsigned int j = 0; j < cols; ++j )
			{
				if ( i == j )
				{
					continue;
				}

				ValueType weight = data( i, j );

				if ( weight > threshold )
				{
					typename ImageType::IndexType index;
					index[0] = i;
					index[1] = j;

					if( !weighted )
						output->SetPixel( index, 1.0 );
					else
						output->SetPixel( index, weight );
				}
			}
		}
	}

	/**
	 * Return maximum value of matrix.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::GetMaxValue( const DataMatrixType& data )
	{
		unsigned int rows = data.rows();
		unsigned int cols = data.cols();

		if ( rows == 0 || cols == 0 )
		{
			std::cerr << "Could not calculate max value: data matrix is empty!" << std::endl;
			exit( EXIT_FAILURE );
		}

		ValueType max = data( 1, 0 );

		for ( unsigned int i = 0; i < rows; ++i )
		{
			for ( unsigned int j = 0; j < cols; ++j )
			{
				if ( i == j )
				{
					continue;
				}

				ValueType weight = data( i, j );

				if ( weight > max )
				{
					max = weight;
				}
			}
		}
		return max;
	}

	/**
	 * Degree from matrix.
	 */
	template< class ValueType >
	void Graph< ValueType >::Degree( VectorType& degrees, const DataMatrixType& data, ValueType threshold, bool totalDegree, bool weighted,
			bool normalize, const std::string& maskImageFileName, const std::string& output )
	{
		unsigned int rows = data.rows();
		unsigned int cols = data.columns();

		degrees.set_size( rows );

		for ( unsigned int i = 0; i < rows; i++ )
		{
			degrees[i] = 0;

			if ( totalDegree )
			{
				for ( unsigned int j = 0; j < cols; j++ )
				{
					ValueType value = data( i, j );
					if ( ( value > threshold ) && ( i != j ) )
					{
						if ( weighted )
							degrees[i] += value;
						else
							degrees[i] += 1.0;
					}
				}
			}
			else // out-degree
			{
				for ( unsigned int j = i + 1; j < cols; j++ )
				{
					ValueType value = data( i, j );
					if ( value > threshold )
					{
						if ( weighted )
							degrees[i] += value;
						else
							degrees[i] += 1.0;
					}
				}
			}
		}

		// ( V - mean ) / sd ...
		if ( normalize )
		{
			ValueType mean = GetMean( degrees );
			ValueType sd = GetSD( degrees );

			for ( unsigned int i = 0; i < degrees.size(); i++ )
			{
				degrees[i] = ( degrees[i] - mean ) / sd;
			}
		}

		// map to mask ...
		if ( !maskImageFileName.empty() )
		{
			std::string out = output + ".nii.gz";
			Fill3DImageWithValues( maskImageFileName, out, degrees );
		}
	}

	/**
	 * Degree,
	 */
	template< class ValueType >
	void Graph< ValueType >::Degree( VectorType& degrees, const std::string& filename, ValueType threshold, bool totalDegree,
			bool weighted, bool normalize, const std::string& maskImageFileName, const std::string& output )
	{

		typename ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( filename );
		reader -> Update();

		typename ImageType::Pointer image = reader -> GetOutput();
		reader = 0;

		ValueType* buffer = image -> GetPixelContainer() -> GetBufferPointer();
		typename ImageType::RegionType region = image -> GetLargestPossibleRegion();
		typename ImageType::SizeType size = region.GetSize();

		unsigned int rows = size[0];
		unsigned int cols = size[1];

		DataMatrixType data( rows, cols, buffer );

		Degree( degrees, data, threshold, totalDegree, weighted, normalize, maskImageFileName, output );
	}

	/**
	 * Return number of image dimensions.
	 */
	template< class ValueType >
	unsigned int Graph< ValueType >::GetImageDimensions( const std::string& inputFileName )
	{
		itk::ImageIOFactory::ImageIOBasePointer io = itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(),
				itk::ImageIOFactory::ReadMode );
		if ( !io )
		{
			std::cerr << "*** ERROR ***: Could not read: " << inputFileName << std::endl;
			exit( EXIT_FAILURE );
		} else
		{
			io -> SetFileName( inputFileName );
			io -> ReadImageInformation();
			return io -> GetNumberOfDimensions();
		}
	}

	/**
	 * Return standard deviation.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::GetSD( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::variance > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return sqrt( variance( acc ) );
	}

	/**
	 * Return mean.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::GetMean( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::mean > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return mean( acc );
	}

	/**
	 * Return median.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::GetMedian( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::median > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return median( acc );
	}

	template< class ValueType >
	ValueType Graph< ValueType >::GetMedian( const std::vector< ValueType >& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::median > > acc;

		for( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return median( acc );
	}


	/**
	 * Return min.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::GetMin( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::min > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return boost::accumulators::min( acc );
	}

	/**
	 * Return max.
	 */
	template< class ValueType >
	ValueType Graph< ValueType >::GetMax( const VectorType& V )
	{
		using namespace boost;
		using namespace accumulators;

		accumulator_set< ValueType, stats< tag::max > > acc;

		for ( unsigned int i = 0; i < V.size(); i++ )
			acc( V[i] );

		return boost::accumulators::max( acc );
	}

	template< class ValueType >
	void Graph< ValueType >::Fill3DImageWithValues( const std::string& referenceImageFileName, const std::string& outputImageFileName,
			const STDVectorType& data )
	{
		VectorType V( data.size() );
		for( unsigned int i = 0; i < data.size(); i++ )
		{
			V[i] = data[i];
		}

		Fill3DImageWithValues( referenceImageFileName, outputImageFileName, V );
	}

	/**
	 * Create new image using given reference image (mask)
	 * and fill all nonzero voxels in the reference image with the given data.
	 */
	template< class ValueType >
	void Graph< ValueType >::Fill3DImageWithValues( const std::string& referenceImageFileName,
			const std::string& outputImageFileName,
			const VectorType& data )
	{
		typedef unsigned char MaskPixelType;
		typedef itk::Image< MaskPixelType, 3 > MaskImageType;
		typedef itk::Image< ValueType, 3 > ImageType;
		typedef itk::ImageMaskSpatialObject< 3 > MaskType;
		typedef itk::ImageSpatialObject< 3, ValueType > SpatialType;
		typedef itk::ImageFileReader< MaskImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( referenceImageFileName.c_str() );
		reader -> Update();

		// mask-image (for dims, etc)
		MaskImageType::Pointer maskImage = reader -> GetOutput();
		MaskType::RegionType region = maskImage -> GetLargestPossibleRegion();
		MaskType::SizeType size = region.GetSize();
		unsigned int voxels = region.GetNumberOfPixels();
		const MaskPixelType* dataMask = maskImage->GetPixelContainer()->GetBufferPointer();

		typename ImageType::Pointer bufferImage = ImageType::New();
		bufferImage->CopyInformation( maskImage );
		bufferImage->SetRegions( region );
		bufferImage->Allocate();
		bufferImage->FillBuffer( 0 );

		ValueType* buffer = bufferImage -> GetPixelContainer() -> GetBufferPointer();

		// disconnect pipe-lines
		reader = 0;

		// fill data buffer with strengths
		unsigned int count = 0;
		for ( unsigned int i = 0; i < voxels; ++i )
		{
			if ( dataMask[i] != 0 )
			{
				buffer[i] = data[count];
				count++;
			} else
			{
				buffer[i] = 0;
			}
		}

		typename WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( outputImageFileName.c_str() );
		writer->SetInput( bufferImage );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputImageFileName << std::endl;
		}
	}

	/**
	 * Create new image using given reference image (mask)
	 * and fill all index voxels in the reference image with the given data, indicated by the indices.
	 */
	template< class ValueType >
	void Graph< ValueType >::Fill3DImageWithValues(
			const std::string& referenceImageFileName,
			const std::string& outputImageFileName,
			const VectorType& indices,
			const VectorType& data )
	{
		if ( indices.size() != data.size() )
		{
			std::cerr << "*** ERROR ***: could not fill 3d image with values (indices.size != data.size)!" << std::endl;
			exit( EXIT_FAILURE );
		}
		else if( referenceImageFileName.empty() | outputImageFileName.empty() )
		{
			std::cerr << "*** ERROR ***: reference or output is empty!" << std::endl;
			exit( EXIT_FAILURE );
		}

		typedef unsigned char MaskPixelType;
		typedef itk::Image< MaskPixelType, 3 > MaskImageType;
		typedef itk::Image< ValueType, 3 > ImageType;
		typedef itk::ImageMaskSpatialObject< 3 > MaskType;
		typedef itk::ImageSpatialObject< 3, ValueType > SpatialType;
		typedef itk::ImageFileReader< MaskImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		ReaderType::Pointer reader = ReaderType::New();
		reader -> SetFileName( referenceImageFileName.c_str() );
		reader -> Update();

		// mask-image (for dims, etc)
		MaskImageType::Pointer maskImage = reader -> GetOutput();
		MaskType::RegionType region = maskImage -> GetLargestPossibleRegion();
		MaskType::SizeType size = region.GetSize();
		unsigned int voxels = region.GetNumberOfPixels();

		const MaskPixelType* dataMask = maskImage->GetPixelContainer()->GetBufferPointer();

		typename ImageType::Pointer bufferImage = ImageType::New();
		bufferImage->CopyInformation( maskImage );
		bufferImage->SetRegions( region );
		bufferImage->Allocate();
		bufferImage->FillBuffer( 0 );

		ValueType* buffer = bufferImage->GetPixelContainer()->GetBufferPointer();

		// disconnect pipe-line
		reader = 0;

		// use map to get data ...
		std::map< unsigned int, ValueType > indicesMap;
		typename std::map< unsigned int, ValueType >::iterator it;

		for( unsigned int i = 0; i< indices.size(); i++ )
		{
			indicesMap[ indices[i] ] = data[i];
		}

		// fill data buffer with strengths
		unsigned int count = 0;
		for ( unsigned int i = 0; i < voxels; ++i )
		{
			if ( dataMask[i] != 0 )
			{
				it = indicesMap.find( count );
				if( it != indicesMap.end() )
					buffer[i] = it->second;

				count++;
			} else
			{
				buffer[i] = 0;
			}
		}

		typename WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( outputImageFileName.c_str() );
		writer->SetInput( bufferImage );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputImageFileName << std::endl;
		}
	}

	/**
	 * Get vector from file.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetVectorFromFile( const std::string& inputFileName, VectorType& timeSeries )
	{

		std::ifstream in( inputFileName.c_str() );
		std::ifstream inLineCount( inputFileName.c_str() );

		if ( !in || !inLineCount )
		{
			std::cerr << "*** ERROR ***: Unable to open file \"" << inputFileName << "\" for time-series input.\n";
			exit( EXIT_FAILURE );
		}

		// if we use vnl_vector< ValueType > as VectorType, there is no push_back,
		// so we allocate the vector first. This means we need to know the
		// number of lines...
		unsigned int numberOfLines = std::count( std::istreambuf_iterator< char >( inLineCount ), std::istreambuf_iterator< char >(), '\n' );

		// allocate ...
		timeSeries.clear();
		timeSeries.set_size( numberOfLines );

		std::string line;
		unsigned int i = 0;

		while ( getline( in, line ) )
		{
			TokType tok( line, boost::char_separator< char >( " " ) );

			// only first value in this case ...
			TokType::iterator id = tok.begin();

			ValueType value = boost::lexical_cast< ValueType >( *id );

			//timeSeries.push_back( value ); <- if VectorType = std::vector< ValueType > ...
			timeSeries[i] = value;
			i++;

			//for ( Tok::iterator id = tok.begin(); id != tok.end(); ++id )
			//{
			//	ValueType value = boost::lexical_cast< ValueType >( *id );
			//}
		}

		in.close();
		inLineCount.close();
	}

	/**
	 * Write vector to filename.
	 */
	template< class ValueType >
	void Graph< ValueType >::WriteVectorToFile( const std::string& outputFileName, const VectorType& V )
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

	/**
	 * Write stdvector to filename.
	 */
	template< class ValueType >
	void Graph< ValueType >::WriteVectorToFile( const std::string& outputFileName, const std::vector< ValueType >& V )
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

	/**
	 * Return name for given key.
	 * And replace all spaces with underscores.
	 * If name is not found, return key as string.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetLabelName( const std::string& root, const std::string& end, const MapType& map, std::string& name, int key )
	{
		MapType::const_iterator it = map.find( key );

		std::stringstream ss;

		// not found ...
		if ( it == map.end() )
		{
			ss << root << "_" << key << end;
			name = ss.str();
		} else
		{
			ss << root << "_" << it->second << end;
			name = ss.str();
			std::replace( name.begin(), name.end(), ' ', '_' );
		}
	}

	/**
	 * Get label names from given file.
	 */
	template< class ValueType >
	void Graph< ValueType >::GetLabelNames( MapType& maps, const std::string& labelNamesFile, const std::string& sep )
	{
		std::ifstream in( labelNamesFile.c_str() );

		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read labels from: " << labelNamesFile << std::endl;
			exit( EXIT_FAILURE );
		}

		std::string line;
		while ( getline( in, line ) )
		{
			TokType tok( line, boost::char_separator< char >( sep.c_str() ) );

			try
			{
				int key = 0;
				std::string value;
				bool firstTime = true;

				for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
				{
					if ( firstTime )
					{
						key = boost::lexical_cast< int >( *id );
						firstTime = false;
					} else
					{
						value += boost::lexical_cast< std::string >( *id );
					}
				}

				maps.insert( std::make_pair( key, value ) );
			} catch ( boost::bad_lexical_cast& e )
			{
				std::cout << "*** WARNING ***: bad lexical cast during label name parsing!" << std::endl;
			}
		}
	}

	// **********************************************

	template< class ValueType >
	void VisabilityGraph< ValueType >::WriteVisabilityMatrix( const VectorType& timeSeries, const std::string& outputFileName,
			bool verbose, bool weighted )
	{
		unsigned int size = timeSeries.size();

		typename ImageType::Pointer outputImage = ImageType::New();
		typename ImageType::RegionType outputRegion;
		typename ImageType::SizeType outputSize;
		outputSize[0] = size;
		outputSize[1] = size;
		outputRegion.SetSize( outputSize );
		outputImage->SetRegions( outputRegion );
		outputImage->Allocate();

		// fill with zeros (no self connectivity)
		outputImage->FillBuffer( itk::NumericTraits< ValueType >::Zero );
		ValueType* outputData = outputImage->GetPixelContainer()->GetBufferPointer();
		DataMatrixType output( size, size, outputData );

		for ( unsigned int i = 0; i < size; i++ )
		{
			for ( unsigned int j = i + 1; j < size; j++ )
			{
				if ( IsVisable( i, j, timeSeries[i], timeSeries[j], timeSeries ) )
				{
					if ( verbose )
					{
						std::cout << "Visible: ";
						std::cout << i << " value: [" << timeSeries[i] << "] -> ";
						std::cout << j << " value: [" << timeSeries[j] << "]";
						std::cout << std::endl;
					}

					ValueType slope = 1;

					// the slope is suggested to be used as 'weight'
					// L. Lacasa, B. Luque, J. Luque, J.C. Nuno; EPL, 86 (2009) 30001
					if ( weighted )
					{
						slope = GetSlope( i, j, timeSeries[i], timeSeries[j] );
					}

					output( i, j ) = slope;
					output( j, i ) = slope;
				}
			}
		}

		typename WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( outputFileName.c_str() );
		writer->SetInput( outputImage );

		try
		{
			writer->Update();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error writing: " << outputFileName << std::endl;
		}
	}

	/**
	 * Write visability matrix from given time-series.
	 */
	template< class ValueType >
	void VisabilityGraph< ValueType >::WriteVisabilityMatrix( const std::string& inputFile, const std::string& outputFileName,
			bool verbose, bool weighted )
	{
		VectorType timeSeries;

		GetVectorFromFile( inputFile, timeSeries );

		WriteVisabilityMatrix( timeSeries, outputFileName, verbose, weighted );
	}

	/**
	 * Return slope.
	 */
	template< class ValueType >
	ValueType VisabilityGraph< ValueType >::GetSlope( unsigned int ta, unsigned int tb, ValueType ya, ValueType yb )
	{
		return std::abs< ValueType >( static_cast< ValueType > ( yb - ya ) / static_cast< ValueType > ( tb - ta ) );
	}

	/**
	 * Return visabilitiy criterion value.
	 */
	template< class ValueType >
	bool VisabilityGraph< ValueType >::IsVisable( unsigned int ta, unsigned int tb, ValueType ya, ValueType yb, const VectorType& vector )
	{
		bool visable = true;

		for ( unsigned int tc = 0; tc < vector.size(); tc++ )
		{
			// only check ta < tc < tb...
			if ( ( ta < tc ) && ( tc < tb ) )
			{
				ValueType yc = vector[tc];

				if ( yc < ( ya + ( yb - ya ) * ( tc - ta ) / ( tb - ta ) ) )
				{
					continue;
				} else
				{
					visable = false;
				}
			}
		}

		return visable;
	}

} // end namespace graph

// TODO is float is enabled, vnl_eigen_system gives trouble -> cast double to ValueType...
//template class graph::Graph< float >;
//template class graph::VisabilityGraph< float >;

template class graph::VisabilityGraph< double >;
template class graph::Graph< double >;

