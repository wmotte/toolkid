#ifndef __graphCommon_h__
#define __graphCommon_h__

#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <map>
#include <utility> // make_pair


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageSpatialObject.h"
#include "itkConvertPixelBuffer.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_diag_matrix.h"
#include "vnl/algo/vnl_generalized_eigensystem.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_sparse_symmetric_eigensystem.h"

#include "flens.h"

#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/accumulators/numeric/functional/complex.hpp>
#include <boost/accumulators/numeric/functional/valarray.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "Kmeans.h"

namespace graph
{
	/**
	 * Graph common.
	 */
	template< class ValueType >
	class Graph
	{

	public:

		// typdefs...
		typedef itk::Image< ValueType, 2 > ImageType;
		typedef itk::Image< ValueType, 3 > Image3DType;
		typedef itk::Image< ValueType, 4 > Image4DType;

		typedef vnl_matrix_ref< ValueType > DataMatrixType;
		typedef typename ImageType::Pointer ImagePointerType;
		typedef typename Image3DType::Pointer Image3DPointerType;
		typedef typename Image4DType::Pointer Image4DPointerType;

		typedef vnl_vector< ValueType > VectorType;
		typedef vnl_matrix< ValueType > MatrixType;

		typedef vnl_diag_matrix< ValueType > DiagMatrixType;

		typedef vnl_vector< double > VectorDoubleType;
		typedef vnl_matrix< double > MatrixDoubleType;
		typedef vnl_diag_matrix< double > DiagMatrixDoubleType;

		typedef std::vector< ValueType > STDVectorType;

		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< ImageType > WriterType;

		typedef boost::tokenizer< boost::char_separator< char > > TokType;
		typedef std::map< int, std::string > MapType;

		typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
		typedef itk::ImageRegionIteratorWithIndex< Image3DType > Iterator3DType;
		typedef itk::ImageRegionIteratorWithIndex< Image4DType > Iterator4DType;

		typedef itk::MinimumMaximumImageFilter< Image4DType > MinMaxCalcType;


		static unsigned int GetMaxValue4DImage( const std::string& inputFileName );
		static void Load3DImage( Image3DPointerType& result, const std::string& filename );
		static void Load4DImage( Image4DPointerType& result, const std::string& filename );
		static void Write3DImage( const Image3DPointerType& output, const std::string& filename );
		static void Write4DImage( const Image4DPointerType& output, const std::string& filename );
		static void GetTimeSeries( MatrixType& matrix, const std::string& inputFileName, const std::string& maskFileName );
		static void GetTimeSeries( MatrixType& matrix, const Image4DPointerType& input, const Image3DPointerType& mask );
		static void Insert3DImageIn4DImage( Image4DPointerType& image4D, const Image3DPointerType& image3D, unsigned int position );
		static void Zero3DImage( Image3DPointerType& output, const std::string& referenceFileName );
		static void Zero4DImage( Image4DPointerType& output, const std::string& referenceFileName, unsigned int size4thDimension );
		static void FillImage( Image3DPointerType& fill, const std::string& referenceFileName, const VectorType& values );
		static void MajorityVote( const std::string& inputFileName, const std::string& maskFileName, const std::string& outputFileName );
		static unsigned int MajorityVote( const VectorType& a, unsigned int nLabels );
		static void Entropy( const std::string& inputFileName, const std::string& maskFileName, const std::string& outputFileName );
		static void Entropy( VectorType& entropy, const MatrixType& timeSeries, unsigned int nLabels );
		static ValueType Entropy( const VectorType& a, unsigned int nLabels );


		static void WriteMatrixToFile( const MatrixType& matrix,
				const std::string& outputFileName );

		/**
		 * KMeans++ clustering.
		 */
		static void KMeansPP( const MatrixType& data, int numberOfClusters, VectorType& clusters );

		/**
		 * Write image to file.
		 */
		static void WriteImageToFile( const ImagePointerType& image,
				const std::string& outputFileName );

		/**
		 * Write eigen-vectors to files.
		 */
		static void WriteEigenFiles( const std::string& outputRoot, unsigned int numberOfEigenVectors,
				const MatrixType& eigenVectors, const VectorType& eigenValues );

		/**
		 * Calculate (normalized) laplacian.
		 */
		static void CalculateLaplacian( MatrixType& L, const MatrixType& G, DiagMatrixType& D, bool normalize );

		/**
		 * Solve eigen-system (dense matrix).
		 */
		static void GeneralizedEigenCalculation(
				const MatrixType& LaplacianMatrix,
				const DiagMatrixType& DegreeMatrix,
				MatrixType& eigenVectors,
				VectorType& eigenValues );

		/**
		 * Solve eigen-system for sparse matrix.
		 */
		static void GeneralizedEigenCalculation(
				vnl_sparse_matrix< ValueType >& LaplacianMatrix,
				MatrixType& eigenVectors,
				VectorType& eigenValues, int numberOfEigenVectors );

		/**
		 * Get max value.
		 */
		static ValueType GetMaxValue( const DataMatrixType& data );

		/**
		 * Get thresholded (binary) matrix.
		 */
		static void GetMatrix( ImagePointerType& output, const DataMatrixType& data,
				ValueType threshold, bool weighted );

		/**
		 * Return matrixtype from image pointer.
		 */
		static void GetMatrix( MatrixType& output, const ImagePointerType& input );

		/**
		 * Get thresholded (binary) matrix, based on normal threshold or fixed k.
		 */
		static void GetMatrix( const ImagePointerType& image, ImagePointerType& output,
				ValueType threshold, ValueType K, bool fixed, bool weighted );

		/**
		 * Read matrix from file.
		 */
		static void ReadMatrix( ImagePointerType& matrix, const std::string& inputFileName );

		/**
		 * Degree.
		 */
		static void Degree( VectorType& degrees, const std::string& filename, ValueType threshold, bool totalDegree, bool weighted,
				bool normalize, const std::string& maskImageFileName, const std::string& output );

		/**
		 * Degree, from matrix type.
		 * maskImageFileName and output are optional.
		 */
		static void Degree( VectorType& degrees, const DataMatrixType& data, ValueType threshold, bool totalDegree, bool weighted,
				bool normalize, const std::string& maskImageFileName, const std::string& output );

		/**
		 * Get vector from file.
		 */
		static void GetVectorFromFile( const std::string& inputFileName, VectorType& V );

		/**
		 * Write vector to file.
		 */
		static void WriteVectorToFile( const std::string& outputFileName, const VectorType& V );
		/**
		 * Write std vector to file.
		 */
		static void WriteVectorToFile( const std::string& outputFileName, const std::vector< ValueType >& V );

		/**
		 * Return image dimensions.
		 */
		static unsigned int GetImageDimensions( const std::string& inputFileName );

		/**
		 * Return standard deviation.
		 */
		static ValueType GetSD( const VectorType& V );

		/**
		 * Return mean.
		 */
		static ValueType GetMean( const VectorType& V );

		/**
		 * Return median.
		 */
		static ValueType GetMedian( const VectorType& V );

		/**
		 * Return median.
		 */
		static ValueType GetMedian( const std::vector< ValueType >& V );

		/**
		 * Return min.
		 */
		static ValueType GetMin( const VectorType& V );

		/**
		 * Return max.
		 */
		static ValueType GetMax( const VectorType& V );

		/**
		 * Fill given reference (mask) image with data in given vector.
		 */
		static void Fill3DImageWithValues( const std::string& referenceImageFileName, const std::string& outputImageFileName,
				const VectorType& data );

		/**
		 * Fill given reference (mask) image with data in given vector.
		 */
		static void Fill3DImageWithValues( const std::string& referenceImageFileName, const std::string& outputImageFileName,
				const STDVectorType& data );

		/**
		 * Fill given reference (mask) image with data as specified by indices.
		 */
		static void Fill3DImageWithValues( const std::string& referenceImageFileName, const std::string& outputImageFileName,
					const VectorType& indices, const VectorType& data );

		/**
		 * Search for label using key.
		 */
		static void GetLabelName( const std::string& root, const std::string& end, const MapType& map, std::string& name, int key );

		/**
		 * Get labels from txt-file.
		 */
		static void GetLabelNames( MapType& maps, const std::string& labelNamesFile, const std::string& sep );
	};

	// ****************************************************************

	/**
	 * Visability Graph extends from Common Graph and converts time-series into
	 * visability matrices.
	 */
	template< class ValueType >
	class VisabilityGraph: public Graph< ValueType >
	{

	public:

		// typdefs...
		typedef	typename Graph< ValueType>::VectorType VectorType;
		typedef typename Graph< ValueType>::ReaderType ReaderType;
		typedef typename Graph< ValueType>::WriterType WriterType;
		typedef typename Graph< ValueType>::ImageType ImageType;
		typedef typename Graph< ValueType>::DataMatrixType DataMatrixType;

		/**
		 * Write visability matrix from given input file.
		 */
		static void WriteVisabilityMatrix( const std::string& inputFile, const std::string& outputFileName, bool verbose, bool weighted );

		/**
		 * Write visability matrix from given time-series vector.
		 */
		static void WriteVisabilityMatrix( const VectorType& inputFile, const std::string& outputFileName, bool verbose, bool weighted );

	protected:

		/**
		 * Return slope.
		 */
		static ValueType GetSlope( unsigned int ta, unsigned int tb, ValueType ya, ValueType yb );

		/**
		 * Check if visable.
		 */
		static bool IsVisable( unsigned int ta, unsigned int tb, ValueType ya, ValueType yb, const VectorType& vector );

	};

} // end namespace graph

#endif /*__graphCommon_h__*/
