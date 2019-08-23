#ifndef __annCommon_h__
#define __annCommon_h__

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
#include "itkAndImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageSpatialObject.h"
#include "itkConvertPixelBuffer.h"
#include "itkPoint.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkDirectedHausdorffDistanceImageFilter.h"


#include <boost/tokenizer.hpp>

#include <boost/lexical_cast.hpp>

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>

#include <boost/filesystem/operations.hpp>

#include <boost/utility.hpp>

#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/accumulators/numeric/functional/complex.hpp>
#include <boost/accumulators/numeric/functional/valarray.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <ANN/ANN.h>

#include "vnl/vnl_matrix.h"

/**
 * ANN Functions.
 *
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC Utrecht, NL.
 * @date: 26-12-2009
 */
namespace ann
{
	/**
	 *
	 */
	template< class ValueType >
	class Ann
	{

	public:

		//typedefs here ...
		static const unsigned int Dimension = 3;
		static const double PI = 3.14159265359;

		typedef unsigned int LabelPixelType;

		typedef itk::Point< ValueType, Dimension > PointType;

		typedef itk::Image< ValueType, Dimension > ImageType;
		typedef itk::Image< LabelPixelType, Dimension > LabelImageType;

		typedef itk::ImageFileReader< ImageType > ImageReaderType;
		typedef itk::ImageFileReader< LabelImageType > LabelImageReaderType;

		typedef itk::ImageFileWriter< ImageType > ImageWriterType;
		typedef itk::ImageFileWriter< LabelImageType > LabelImageWriterType;

		typedef itk::AndImageFilter< LabelImageType, LabelImageType, LabelImageType > LabelAndImageFilterType;
		typedef itk::OrImageFilter< LabelImageType, LabelImageType, LabelImageType > LabelOrImageFilterType;

		typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ImageIteratorType;
		typedef itk::ImageRegionConstIteratorWithIndex< LabelImageType > LabelImageIteratorType;

		typedef std::vector< ValueType > VectorType;
		typedef std::vector< VectorType > MatrixType;

		typedef	typename ImageType::Pointer ImagePointerType;
		typedef typename LabelImageType::Pointer LabelImagePointerType;

		typedef itk::BinaryThresholdImageFilter< ImageType, LabelImageType > BinaryThresholdImageFilterType;
		typedef itk::ThresholdImageFilter< ImageType > ThresholdImageFilterType;
		typedef itk::ConnectedComponentImageFilter< LabelImageType, LabelImageType > ConnectedComponentImageFilterType;
		typedef itk::MaskImageFilter< ImageType, LabelImageType, ImageType > MaskImageFilterType;
		typedef itk::LabelStatisticsImageFilter< ImageType, LabelImageType > LabelStatisticsImageFilterType;

		typedef itk::RelabelComponentImageFilter< LabelImageType, LabelImageType > RelabelComponentImageFilterType;
		typedef itk::SubtractImageFilter< LabelImageType, LabelImageType, LabelImageType > SubtractImageFilterType;
		typedef itk::CastImageFilter< LabelImageType, ImageType > CastImageFilterType;
		typedef itk::SignedMaurerDistanceMapImageFilter< ImageType, ImageType > SignedMaurerDistanceMapImageFilterType;
		typedef itk::DirectedHausdorffDistanceImageFilter< LabelImageType, LabelImageType > DirectedHausdorffDistanceImageFilterType;
		typedef itk::InvertIntensityImageFilter< LabelImageType, LabelImageType > InvertIntensityLabelImageFilterType;

		typedef boost::tokenizer< boost::char_separator< char> > TokType;
		typedef boost::mt19937 random_number_type;

		typedef boost::uniform_int< long> int_distribution_type;
		typedef boost::variate_generator< random_number_type&, int_distribution_type> int_generator_type;

		/**
		 * Convert vnl matrix to std matrix.
		 */
		static void ConvertVNL2STD( const vnl_matrix< ValueType >& G, MatrixType& stdG );

		/**
		 * Convert vnl vector to std vector.
		 */
		static void ConvertVNL2STD( const vnl_vector< ValueType >& G, VectorType& stdG );

		/**
		 * Load vector from txt file.
		 */
		static void LoadVector( VectorType& vector, const std::string& inputFileName );

		/**
		 * Load matrix from 2D image.
		 */
		static void LoadMatrix( MatrixType& result, const std::string& matrixFileName );

		/**
		 * Modularity.
		 */
		static ValueType Modularity( const VectorType& clusters, const MatrixType& matrix );

		/**
		 * Invert label image.
		 */
		static void Invert( const LabelImagePointerType& input, LabelImagePointerType& output );

		/**
		 * Get image dimensions.
		 */
		static unsigned int GetImageDimensions( const std::string& imageFileName );

		/**
		 * Subtract.
		 */
		static void Subtract( LabelImagePointerType& result,
				const LabelImagePointerType& image1,
				const LabelImagePointerType& image2,
				const std::string& outputFileName );
		/**
		 * Distance map.
		 */
		static void DistanceMap(
				ImagePointerType& distanceMap,
				const LabelImagePointerType& ref,
				const std::string& outputFileName );

		/**
		 * Average value.
		 */
		static ValueType AverageValue( const LabelImagePointerType& labels,
				const ImagePointerType& values );

		/**
		 * Return thresholded image.
		 */
		static void ThresholdToBinaryImage( const ImagePointerType& input, ValueType threshold, LabelImagePointerType& output );

		/**
		 * Return thresholded image.
		 */
		static void ThresholdToImage( const ImagePointerType& input, ValueType threshold, ImagePointerType& output );

		/**
		 * Remove isolated components with voxel less than minimal voxel number specified.
		 */
		static void RemoveIsolatedComponents( const LabelImagePointerType& input, LabelImagePointerType& output, unsigned int minVoxels );

		/**
		 * Remove isolated components with voxel less than minimal voxel number specified.
		 */
		static void RemoveIsolatedComponents( const ImagePointerType& input, ImagePointerType& output, unsigned int minVoxels );


		/**
		 * AND operator binary images.
		 */
		static void LabelAnd( LabelImageType::Pointer& result, const LabelImageType::Pointer image1, const LabelImageType::Pointer image2 );

		/**
		 * Return image pointer.
		 */
		static void GetImage( ImagePointerType& result, const std::string& imageFileName );

		/**
		 * Remove isolated components.
		 */
		static void RemoveIsolatedComponents( const std::string& inputFileName, const std::string& outputFileName, unsigned int minVoxels );

		/**
		 * Hausdorff distance range.
		 */
		static void HausdorffDistanceRange( VectorType& average, VectorType& directed, const LabelImagePointerType& ref,
				const ImagePointerType& seg, unsigned int minVoxels );


		static void HausdorffDistance(
				ValueType& negative,
				ValueType& positive,
				const LabelImagePointerType& ref,
				const LabelImagePointerType& seg );

		/**
		 * Haussdorf Distance.
		 */
		static void HausdorffDistance( ValueType& negative, ValueType& positive,
				const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels );

		/**
		 * Range of border distances.
		 */
		static void BorderDistanceRange( VectorType& negative, VectorType& positive, const LabelImagePointerType& ref, const ImagePointerType& seg, unsigned int minVoxels );

		/**
		 * Range of border distances.
		 */
		static void BorderDistanceRange( VectorType& negative, VectorType& postive, const std::string& ref, const std::string& seg, unsigned int minVoxels );

		/**
		 * Hausdorff distance range.
		 */
		static void HausdorffDistanceRange( VectorType& average, VectorType& directed, const std::string& ref, const std::string& seg, unsigned int minVoxels );

		/**
		 * Border distance (positive and negative).
		 */
		static void BorderDistance(
				ValueType& positive,
				ValueType& negative,
				const LabelImagePointerType& ref,
				const LabelImagePointerType& seg,
				const std::string& distanceMapFileName,
				const std::string& distanceMapInverseFileName );

		/**
		 * Border distance (positive and negative).
		 */
		static void BorderDistance(
				ValueType& positive,
				ValueType& negative,
				const std::string& refFileName,
				const std::string& segFileName, unsigned int minVoxels,
				const std::string& distanceMapFileName,
				const std::string& distanceMapInverseFileName );

		/**
		 * SI.
		 */
		static ValueType SI( const std::string& refFileName, const std::string& segFileName );

		/**
		 * PSI.
		 */
		static ValueType PSI( const std::string& refFileName, const std::string& segFileName );

		/**
		 * Sensitivity.
		 */
		static ValueType Sensitivity( const std::string& refFileName, const std::string& segFileName );

		/**
		 * Specificity.
		 */

		static ValueType Specificity( const std::string& refFileName,
				const std::string& segFileName, const std::string& maskFileName );

		/**
		 * Read images into vector of image pointers.
		 */
		static void ReadImages( std::vector< ImagePointerType>& images,
				const std::vector< std::string>& inputFileNames );

		/**
		 * Read label image from file.
		 */
		static void ReadImage( LabelImagePointerType& image, const std::string& imageFileName );

		/**
		 * Read image from file.
		 */
		static void ReadImage( ImagePointerType& image, const std::string& imageFileName );

		/**
		 * Add spatial features (x,y,z positions relative to COG) to input matrix.
		 * Mask is needed to calculate COG.
		 */
		static void FillSpatialData( MatrixType& features,
				const std::string& labelFileName,
				const std::string& maskFileName );

		/**
		 * Add euclidean distance to input matrix.
		 */
		static void FillEuclideanData( MatrixType& features,
				const std::string& labelFileName, const std::string& maskFileName );

		/**
		 * Add angles to input matrix.
		 */
		static void FillAnglesData( MatrixType& features,
				const std::string& labelFileName, const std::string& maskFileName );

		/**
		 * Write matrix to file.
		 */
		static void WriteData( const MatrixType& results, const std::string& outputFileName );

		/**
		 * Write vector to file.
		 */
		static void WriteData( const VectorType& results, const std::string& outputFileName );

		/**
		 * Write label image to file.
		 */
		static void WriteData( const LabelImagePointerType& output, const std::string& outputFileName );

		/**
		 * Write image to file.
		 */
		static void WriteData( const ImagePointerType& output, const std::string& outputFileName );

		/**
		 * Add intensity features to input matrix.
		 */
		static void FillIntensityData(
				MatrixType& features,
				const std::vector< ImagePointerType>& images,
				const LabelImagePointerType& labels );

		/**
		 * Add intensity features to input matrix (inputs are files).
		 */
		static void FillIntensityData(
				MatrixType& features,
				const std::vector< std::string>& images,
				const std::string& labels );

		/**
		 * Return vector with class labels.
		 */
		static void FillClassLabels(
				MatrixType& features,
				const LabelImagePointerType& labels );

		/**
		 * Normalize vector.
		 */
		static void Normalize( VectorType& V );

		/**
		 * Normalize matrix.
		 */
		static void Normalize( MatrixType& V );

		/**
		 * Project matrix vectors to images.
		 */
		static void Project2Image( const MatrixType& M, const std::string& maskFileName,
				const std::string outputFileName, const std::vector< int>& ids );

		/**
		 * Project vector to image.
		 */
		static void Project2Image( const VectorType& V, const std::string& maskFileName,
				const std::string& outputFileName );

		/**
		 * Read training data from file.
		 */
		static void ReadTrainingData( const std::string& trainingDataFileName,
				MatrixType& trainingData, VectorType& trainingClasses, double percentage );

		/**
		 * Query data using training data.
		 */
		static void Query( MatrixType& trainingData, VectorType& trainingClasses, MatrixType& queryData, MatrixType& queryResultMatrix,
				unsigned int k, double errorBound, std::vector< int > probabilistic, const std::string& dumpFileName, bool prioritySearch );

		/**
		 * Return ROC AUC (sensitivity as first vector, 1-specificity as second vector, ROC as third vector).
		 */
		static ValueType ROC( VectorType& sensitivity, VectorType& specificity, VectorType& roc, const std::string& ref, const std::string& seg, const std::string& mask, unsigned int minVoxels );

		/**
		 * Range of SI.
		 */
		static void RangeSI( VectorType& si, const std::string& refFileName, const std::string& segFileName, unsigned int minVoxels );

		/**
		 * Return image pointer.
		 */
		static void GetLabelImage( LabelImagePointerType& result, const std::string& imageFileName );


	protected:

		/**
		 * Return angle in degrees (-180 <-> 180).
		 */
		static ValueType Angle( ValueType x, ValueType y );

		/**
		 * SI.
		 */
		static ValueType SI( const LabelImagePointerType& ref, const LabelImagePointerType& seg );

		/**
		 * PSI.
		 */
		static ValueType PSI( const LabelImagePointerType& ref, const ImagePointerType& seg );

		/**
		 * Return AUC.
		 */
		static ValueType AUC( const VectorType& sensitivity, const VectorType& specificity, VectorType& roc );

		/**
		 * True positives.
		 */
		static ValueType TP( const LabelImagePointerType& ref, const LabelImagePointerType& seg );

		/**
		 * True negatives.
		 */
		static ValueType TN( const LabelImagePointerType& ref, const LabelImagePointerType& seg, const LabelImagePointerType& mask );

		/**
		 * False negatives.
		 */
		static ValueType FN( const LabelImagePointerType& ref, const LabelImagePointerType& seg );

		/**
		 * False positives.
		 */
		static ValueType FP( const LabelImagePointerType& ref, const LabelImagePointerType& seg );

		/**
		 * Return number of nonzero voxels.
		 */
		static unsigned int GetNonZeroVoxels( const LabelImageType::Pointer image );

		/**
		 * Return number + value of nonzero voxels.
		 */
		static ValueType GetNonZeroVoxels( const ImagePointerType image );

		/**
		 * OR operator binary images.
		 */
		static void LabelOr( LabelImageType::Pointer& result, const LabelImageType::Pointer image1, const LabelImageType::Pointer image2 );

		/**
		 * Return AND voxels.
		 */
		static ValueType And( const LabelImagePointerType ref, const ImagePointerType seg );

		/**
		 * Check Dimensions.
		 */
		static void CheckDimensions( const std::string& image );

		/**
		 * Check Dimensions.
		 */
		static void CheckDimensions( const std::vector< std::string>& images );

		/**
		 * Return center-of-gravity.
		 */
		static void GetCog( PointType& point, const std::string& imageFileName );

		/**
		 * Return euclidean distance for given points.
		 */
		static ValueType EuclideanDistance( const PointType& A, const PointType& B );

		/**
		 * Return sd for given vector.
		 */
		static ValueType GetSD( VectorType& V );

		/**
		 * Resize training data matrix.
		 */
		static void ResizeTrainingData( MatrixType& trainingData, VectorType& trainingClasses, double percentage );

		/**
		 * Return vector with total random values, taken from inputs.
		 */
		static void GetRandomValues( const VectorType& inputs, random_number_type& generator, VectorType& results,
				unsigned int total );

		/**
		 * Return number of training dimensions.
		 */
		static unsigned int GetNumberOfTrainingDims( const std::string& trainingData );

		/**
		 * Return number of training points.
		 */
		static unsigned int GetNumberOfTrainingPoints( const std::string& trainingData );

		/**
		 * Dump tree to file (DEBUG).
		 */
		static void DumpTree( const std::string& dumpFileName, ANNkd_tree* kdTree );

		/**
		 * Get probability of given ANN-point.
		 */
		static ValueType Probabilistic( const VectorType& trainingClasses, unsigned int classNumber, unsigned int k, const ANNidxArray& nnIdx );

		/**
		 * Get KD-tree from given matrix.
		 */
		static ANNkd_tree* GetKDTree( const MatrixType& trainingData );

		/**
		 * Specificity.
		 */
		static ValueType Specificity( const LabelImagePointerType& ref, const LabelImagePointerType& seg,
				const LabelImagePointerType& mask );

		/**
		 * Sensitivity.
		 */
		static ValueType Sensitivity( const LabelImagePointerType& ref, const LabelImagePointerType& seg );
	};

} // end namespace ann

#endif /*__annCommon_h__*/
