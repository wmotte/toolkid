#ifndef __dtifit_h___
#define __dtifit_h___

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include "vnl/algo/vnl_qr.h"

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkAutoPointer.h"
#include "itkRGBPixel.h"
#include "itkDiffusionTensor3D.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

#include "procparser.h"


class Procparser;

namespace dtifit
{
	class AlgoEnum
	{
	public:
		enum Enum
		{
			WeightedLinearLeastSquares, ConstrainedWeightedLinearLeastSquares, ConstrainedWeightedNonLinearLeastSquares
		};
	private:
		Enum mVal;
	public:
		AlgoEnum();
		AlgoEnum( const Enum& e )
		{
			mVal = e;
		}
		bool operator==( const Enum& e )
		{
			if( e == mVal )
				return true;
			else
				return false;
		}
		AlgoEnum& operator=( const Enum& e );
	};

	/**
	 * This class implements several fitting routines for DTI data, including the
	 * constrained nonlinear least squares DTI fit, using Newton's method.
	 *
	 * Reference:
	 * (1) Koay CG, Chang LC, Carew JD, Pierpaoli C, Basser PJ.
	 * 		A unifying theoretical and algorithmic framework for least squares methods of
	 * 		estimation in diffusion tensor imaging. Journal of Magnetic Resonance.
	 * 		182 (2006) 115-125.
	 */
	class DTIFit: public itk::LightObject
	{
	public:

		typedef double PixelType;
		typedef itk::Image< PixelType, 4 > InputImageType;


		typedef vnl_vector< PixelType > VectorType;
		typedef vnl_matrix< PixelType > MatrixType;
		typedef itk::Image< PixelType, 3 > OutputImageType;
		typedef itk::Image< PixelType, 4 > OutputVectorImageType;

		typedef itk::ImageFileReader< OutputImageType > ImageReaderType;

		typedef itk::DiffusionTensor3D< double > TensorPixelType;
		typedef itk::Image< TensorPixelType, 3 > TensorImageType;

		typedef TensorPixelType::EigenValuesArrayType EigenValuesArrayType;
		typedef TensorPixelType::EigenVectorsMatrixType EigenVectorsMatrixType;

		typedef itk::RGBPixel< unsigned char > RGBPixelType;
		typedef itk::Image< RGBPixelType, 3 > RGBImageType;
		typedef vnl_matrix< PixelType > DesignMatrixType;

		typedef itk::ImageLinearConstIteratorWithIndex< InputImageType > ConstIterator4DType;
		typedef itk::ImageRegionConstIteratorWithIndex< OutputImageType > ConstIterator3DType;
		typedef itk::ImageRegionIteratorWithIndex< TensorImageType > IteratorTensorType;

		typedef itk::ImageFileWriter< OutputImageType > WriterType;
		typedef itk::ImageFileWriter< OutputVectorImageType > VectorWriterType;
		typedef itk::ImageFileWriter< TensorImageType > TensorWriterType;

		typedef DTIFit Self;
		typedef itk::LightObject Superclass;
		typedef itk::SmartPointer< Self > Pointer;
		typedef itk::SmartPointer< const Self > ConstPointer;

		itkNewMacro( Self );
		itkTypeMacro( DTIFit, LightObject );

		/**
		 * Two representation are used:
		 * (1) The ordinary tensor representation (Gamma).
		 * (2) Cholesky representation (Rho).
		 */
		MatrixType ModifiedCholesky( const MatrixType& A );
		VectorType GetCholeskyComposition( const MatrixType& A );

		/* Weighted linear */
		VectorType WLLS( const DTIFit::DesignMatrixType& W, const VectorType si );
		/* Constrained weighted linear */
		VectorType CWLLS( const MatrixType& W, const VectorType& si, bool gammaOutput );
		VectorType CWLLS( const VectorType& gamma, bool gammaOutput );
		/* Constrained Weighted nonlinear */
		VectorType CWNLLS( const DTIFit::DesignMatrixType& W, const VectorType si );

		/* Newton update methods */
		PixelType Update( unsigned int dim, unsigned int N, const MatrixType& W, const VectorType& si, const VectorType& rho, PixelType lambda, VectorType& negGrad, MatrixType& hessian );
		PixelType ValueAt( const MatrixType& W, const VectorType& si, unsigned int dim, unsigned int N, const VectorType& rho );

		/* Jacobian/Hessian methods */
		MatrixType J( const VectorType& rho, unsigned int dim );
		MatrixType WJ( const VectorType& rho, const MatrixType& W, unsigned int dim );
		MatrixType Hessian( const VectorType& sihat, const VectorType&  r, unsigned int dim, unsigned int N, const MatrixType& W );

		/* Matrix and vector manipulations */
		VectorType Abs( const VectorType& V );
		MatrixType Abs( const MatrixType& A );
		MatrixType RemoveDiagonals( const MatrixType& A );
		VectorType GetDiagonals( const MatrixType& mat );
		MatrixType Dot( const MatrixType& a, const MatrixType& b );
		PixelType Dot( const VectorType& a, const VectorType& b );
		VectorType Multiply( const VectorType& a, const VectorType& b );
		PixelType Max( PixelType a, PixelType b );
		MatrixType GetIdentity( unsigned int n );

		void Gamma( const VectorType& rho, VectorType& gamma );
		bool IsSymmetric( const MatrixType& mat );
		double MachineEpsilon();

		void SetAlgorithm( AlgoEnum::Enum algorithm );
		void SetImage( InputImageType::Pointer input );
		void SetMask( OutputImageType::Pointer mask );
		void SetDesignMatrix( const DesignMatrixType& designMatrix );
		DesignMatrixType GetDesignMatrixVersion1( const MatrixType& g, const VectorType& B );
		DesignMatrixType GetDesignMatrixVersion2( const MatrixType& g, const VectorType& B );

		void Fit();
		VectorType FitVoxel( const VectorType& si );
		TensorPixelType Tensor( const VectorType& gamma );

		OutputImageType::Pointer GetM( const std::string& maskFileName, InputImageType::ConstPointer image ); // const std::string& maskFileName, ImageType::ConstPointer img );
		TensorImageType::Pointer EmptyTensorImage( OutputImageType::ConstPointer input );

		OutputImageType::Pointer GetFA();
		OutputImageType::Pointer GetRA();
		OutputImageType::Pointer GetTrace();

		OutputImageType::Pointer GetL1();
		OutputImageType::Pointer GetL2();
		OutputImageType::Pointer GetL3();

		OutputVectorImageType::Pointer GetV1();
		OutputVectorImageType::Pointer GetV2();
		OutputVectorImageType::Pointer GetV3();

		TensorImageType::Pointer GetTensor();
		RGBImageType::Pointer GetColorFA();

		void Run( InputImageType::Pointer input, const std::string& maskFileName, const Procparser& pp, const std::string& outputBase, bool mirrorGradientTable = false, unsigned int algorithm = 3 );

	protected:
		DTIFit();
		~DTIFit();

	private:

		struct Impl;
		itk::AutoPointer< Impl > m_Impl;

		DTIFit( const Self& );
		void operator=( const Self& );
	};

} // end namespace dtifit

#endif /*__dtifit_h___*/
