#ifndef __dkifitCommon_h__
#define __dkifitCommon_h__

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "vcl_algorithm.h"
#include "vcl_functional.h"

#include <cmath>
#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace dki
{
	typedef double PixelType;
	typedef unsigned char MaskPixelType;

	typedef itk::Image< PixelType, 4 > ImageType;
	typedef itk::Image< MaskPixelType, 3 > MaskImageType;
	typedef itk::Image< PixelType, 4 > OutputImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileReader< MaskImageType > MaskReaderType;

	typedef itk::ImageLinearConstIteratorWithIndex< ImageType > ConstIterator4DType;
	typedef itk::ImageLinearIteratorWithIndex< ImageType > Iterator4DType;
	typedef itk::ImageRegionIteratorWithIndex< MaskImageType > Iterator3DType;
	typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType > ConstIterator3DType;

	typedef itk::ImageFileWriter< OutputImageType > WriterType;

	typedef boost::tokenizer< boost::char_separator< char > > TokType;

	typedef vnl_vector< PixelType > VectorType;
	typedef vnl_matrix< PixelType > MatrixType;
	typedef vnl_diag_matrix< PixelType > DiagMatrixType;

	typedef vnl_vector< unsigned int > IndexVectorType;

	/**
	 * Arguments.
	 */
	struct parameters
	{
		std::string inputFileName;
		std::string maskFileName;
		std::string bvecsFileName;
		std::string bvalsFileName;
		std::string outputFileName;
		std::string sep;
	};

	/**
	 * Base class for diffusion kurtosis fitting routines.
	 *
	 * Prepares: design matrix, etc.
	 * wim@invivonmr.uu.nl; 30-03-2011
	 */
	class DkiFit
	{
	public:

		// how many output images (6 dti elements + 15 dki elements).
		const static unsigned int m_OUTPUT_IMAGES = 21;

		DkiFit( const parameters& args );

	protected:

		ImageType::Pointer m_Input;
		MaskImageType::Pointer m_Mask;

		MatrixType m_A_D;
		MatrixType m_A_K;
		MatrixType m_DesignMatrix;

		unsigned int m_M; // number of non-zero dwi images
		VectorType m_bvals; // non-zero bvals
		MatrixType m_bvecs; // non-zero bvecs
		IndexVectorType m_nonzero_indices; // non-zero indices
		OutputImageType::Pointer m_Output;
		std::string m_sep;

		void InitData( const std::string& inputFileName );
		void InitMask( const std::string& maskFileName );
		void Write( const std::string& outputFileName );
		void ReadBVals( const std::string& bvalsFileName );
		void ReadBVecs( const std::string& bvecsFileName );
		unsigned int GetNumberOfRows( const std::string& inputFile );
		void InitA_D();
		void InitA_K();
		void InitDesignMatrix();
		void RemainDWIOnly();
		VectorType LogDWIOnlyData( const VectorType data );
		PixelType Median( const VectorType& V );
		MatrixType CheckDeterminant( const MatrixType& XWX );
		VectorType Residuals( const MatrixType& A, const DiagMatrixType& W, const VectorType& b, const VectorType& x );
		VectorType Residuals( const MatrixType& A, const VectorType& b, const VectorType& x );

		DkiFit(); // prevent construction without input arguments.
	};

} // end namespace dki

#endif /*__dkifitCommon_h__*/
