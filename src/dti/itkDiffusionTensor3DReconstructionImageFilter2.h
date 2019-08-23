/*=========================================================================

 Changed by wim@invivonmr.uu.nl:

 =========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionImageFilter2_h
#define __itkDiffusionTensor3DReconstructionImageFilter2_h

#include "itkImageToImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/algo/vnl_svd.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"

namespace itk
{
	/** \class DiffusionTensor3DReconstructionImageFilter
	 * \brief This class takes as input one or more reference image (acquired in the 
	 * absence of diffusion sensitizing gradients) and 'n' diffusion
	 * weighted images and their gradient directions and computes an image of 
	 * tensors. (with DiffusionTensor3D as the pixel type). Once that is done, you 
	 * can apply filters on this tensor image to compute FA, ADC, RGB weighted 
	 * maps etc. 
	 *
	 * \par Inputs and Usage
	 * There are two ways to use this class. When you have one reference image and \c n
	 * gradient images, you would use the class as
	 * \code
	 *       filter->SetReferenceImage( image0 );
	 *       filter->AddGradientImage( direction1, image1 );
	 *       filter->AddGradientImage( direction2, image2 );
	 *   ...
	 * \endcode
	 *
	 * \par
	 * When you have the 'n' gradient and one or more reference images in a single 
	 * multi-component image (VectorImage), you can specify the images simply as
	 * \code
	 *       filter->SetGradientImage( directionsContainer, vectorImage );
	 * \endcode
	 * Note that this method is used to specify both the reference and gradient images.
	 * This is convenient when the DWI images are read in using the 
	 * <a href="http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format">NRRD</a> 
	 * format. Like the Nrrd format, the reference images are those components of the 
	 * vectorImage whose gradient direction is (0,0,0). If more than one reference image
	 * is present, they are averaged prior to applying the Stejskal-Tanner equations.
	 *
	 * \par Outputs
	 * The output image is an image of Tensors:
	 * \code
	 *       Image< DiffusionTensor3D< TTensorPixelType >, 3 >
	 * \endcode
	 *
	 * \par Parameters
	 * \li Threshold -  Threshold on the reference image data. The output tensor will 
	 * be a null tensor for pixels in the reference image that have a value less 
	 * than this.
	 * \li BValue - See the documentation of SetBValue().
	 * \li At least 6 gradient images must be specified for the filter to be able 
	 * to run.
	 * 
	 * 
	 * \par Template parameters
	 * The class is templated over the pixel type of the reference and gradient 
	 * images (expected to be scalar data types) and the internal representation
	 * of the DiffusionTensor3D pixel (double, float etc).
	 *  
	 * \par References:
	 * \li<a href="http://lmi.bwh.harvard.edu/papers/pdfs/2002/westinMEDIA02.pdf">[1]</a> 
	 * <em>C.F.Westin, S.E.Maier, H.Mamata, A.Nabavi, F.A.Jolesz, R.Kikinis,
	 * "Processing and visualization for Diffusion tensor MRI", Medical image
	 * Analysis, 2002, pp 93-108.</em>
	 * \li<a href="splweb.bwh.harvard.edu:8000/pages/papers/westin/ISMRM2002.pdf">[2]</a>
	 * <em>A Dual Tensor Basis Solution to the Stejskal-Tanner Equations for DT-MRI</em>
	 * 
	 * \par WARNING:
	 * Although this filter has been written to support multiple threads, please 
	 * set the number of threads to 1.
	 * \code
	 *         filter->SetNumberOfThreads(1);
	 * \endcode
	 * This is due to buggy code in netlib/dsvdc, that is called by vnl_svd. 
	 * (used to compute the psudo-inverse to find the dual tensor basis).
	 *
	 * \author Thanks to Xiaodong Tao, GE, for contributing parts of this class. Also
	 * thanks to Casey Goodlet, UNC for patches to support multiple baseline images
	 * and other improvements.
	 * 
	 * \note
	 * This work is part of the National Alliance for Medical image Computing 
	 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
	 * for Medical Research, Grant U54 EB005149.
	 *
	 * \par Examples and Datasets
	 * See Examples/Filtering/DiffusionTensor3DReconstructionImageFilter.cxx
	 * Sample DTI datasets may be obtained from 
	 \begin verbatim
	 ftp://public.kitware.com/pub/namic/DTI/Data/dwi.nhdr
	 ftp://public.kitware.com/pub/namic/DTI/Data/dwi.img.gz ( gunzip this )
	 \end verbatim
	 *
	 * \sa DiffusionTensor3D SymmetricSecondRankTensor 
	 * \ingroup Multithreaded  TensorObjects
	 */

	template< class TReferenceImagePixelType, class TGradientImagePixelType = TReferenceImagePixelType, class TTensorPixelType = double >
	class ITK_EXPORT DiffusionTensor3DReconstructionImageFilter: public ImageToImageFilter< Image< TReferenceImagePixelType, 3 > , Image<
			DiffusionTensor3D< TTensorPixelType > , 3 > >
	{

	public:

		typedef DiffusionTensor3DReconstructionImageFilter Self;
		typedef SmartPointer< Self > Pointer;
		typedef SmartPointer< const Self > ConstPointer;
		typedef ImageToImageFilter< Image< TReferenceImagePixelType, 3 > , Image< DiffusionTensor3D< TTensorPixelType > , 3 > > Superclass;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Runtime information support. */
		itkTypeMacro(DiffusionTensor3DReconstructionImageFilter, ImageToImageFilter);

		typedef TReferenceImagePixelType ReferencePixelType;

		typedef TGradientImagePixelType GradientPixelType;

		typedef DiffusionTensor3D< TTensorPixelType > TensorPixelType;

		/** Reference image data,  This image is aquired in the absence 
		 * of a diffusion sensitizing field gradient */
		typedef typename Superclass::InputImageType ReferenceImageType;

		typedef Image< TensorPixelType, 3 > TensorImageType;

		typedef TensorImageType OutputImageType;

		typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

		/** Holds design matrix columns */
		typedef vnl_matrix< double > DesignMatrixType;

		/** Another set method to add a gradient directions and its corresponding
		 * image. The image here is a VectorImage. The user is expected to pass the 
		 * gradient directions in a container. The ith element of the container 
		 * corresponds to the gradient direction of the ith component image the 
		 * VectorImage.  For the baseline image, a vector of all zeros
		 * should be set. */
		void SetGradientImage( const GradientImagesType *image );

		/**
		 * Design matrix set.
		 */
		void SetDesignMatrix( const DesignMatrixType & );

		/** Set method to set the reference image. */
		void SetReferenceImage( ReferenceImageType *referenceImage )
		{

			this->ProcessObject::SetNthInput( 0, referenceImage );
		}

		/** Get reference image */
		virtual ReferenceImageType * GetReferenceImage()
		{
			return ( static_cast< ReferenceImageType * > ( this->ProcessObject::GetInput( 0 ) ) );
		}

		itkSetMacro( Algorithm, unsigned int );

#ifdef ITK_USE_CONCEPT_CHECKING
		/** Begin concept checking */
		itkConceptMacro(ReferenceEqualityComparableCheck,
				(Concept::EqualityComparable<ReferencePixelType>));
		itkConceptMacro(TensorEqualityComparableCheck,
				(Concept::EqualityComparable<TensorPixelType>));
		itkConceptMacro(GradientConvertibleToDoubleCheck,
				(Concept::Convertible<GradientPixelType, double>));
		itkConceptMacro(DoubleConvertibleToTensorCheck,
				(Concept::Convertible<double, TensorPixelType>));
		itkConceptMacro(GradientReferenceAdditiveOperatorsCheck,
				(Concept::AdditiveOperators<GradientPixelType, GradientPixelType,
						ReferencePixelType>));
		itkConceptMacro(ReferenceOStreamWritableCheck,
				(Concept::OStreamWritable<ReferencePixelType>));
		itkConceptMacro(TensorOStreamWritableCheck,
				(Concept::OStreamWritable<TensorPixelType>));
		/** End concept checking */
#endif

	protected:

		DiffusionTensor3DReconstructionImageFilter();
		~DiffusionTensor3DReconstructionImageFilter()
		{};
		void PrintSelf( std::ostream& os, Indent indent ) const;
		void BeforeThreadedGenerateData();
		void ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread, int );

	private:

		/** container to hold design matrix */
		DesignMatrixType m_DesignMatrixContainer;

		/** Number of baseline images */
		unsigned int m_NumberOfBaselineImages;

		/** Algorithm to use for fitting (.e.g linear, weighted linear, nonlinear, constained nonlinear */
		unsigned int m_Algorithm;
	};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionTensor3DReconstructionImageFilter2.txx"
#endif

#endif
