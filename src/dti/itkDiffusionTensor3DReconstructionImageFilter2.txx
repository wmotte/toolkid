/*=========================================================================

 Changed by wim@invivonmr.uu.nl:

 =========================================================================*/
#ifndef __itkDiffusionTensor3DReconstructionImageFilter2_txx
#define __itkDiffusionTensor3DReconstructionImageFilter2_txx

#include "itkDiffusionTensor3DReconstructionImageFilter2.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "vnl/vnl_vector.h"

namespace itk
{

	/**
	 * Constructor.
	 */
	template< class TReferenceImagePixelType, class TGradientImagePixelType, class TTensorPixelType >
	DiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType, TGradientImagePixelType, TTensorPixelType >
	::DiffusionTensor3DReconstructionImageFilter()
	{
		// At least 1 inputs is necessary for a vector image For images added one at a time we need at least six
		this->SetNumberOfRequiredInputs( 1 );
		m_Algorithm = 1;
	}

	/**
	 * Threaded generate data (warning, potential bug).
	 */
	template< class TReferenceImagePixelType, class TGradientImagePixelType, class TTensorPixelType >
	void DiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType, TGradientImagePixelType, TTensorPixelType >
	::ThreadedGenerateData(
			const OutputImageRegionType& outputRegionForThread, int )
	{
		typename OutputImageType::Pointer outputImage = static_cast< OutputImageType * > ( this->ProcessObject::GetOutput( 0 ) );
		typedef ImageRegionConstIterator< GradientImagesType > GradientIteratorType;
		typedef typename GradientImagesType::PixelType GradientVectorType;
		typename GradientImagesType::Pointer gradientImagePointer = NULL;

		ImageRegionIterator< OutputImageType > oit( outputImage, outputRegionForThread );
		oit.GoToBegin();

		gradientImagePointer = static_cast< GradientImagesType * > ( this->ProcessObject::GetInput( 0 ) );

		GradientIteratorType git( gradientImagePointer, outputRegionForThread );
		git.GoToBegin();

		std::vector< unsigned int > baselineind; // contains the indices of the baseline images.

		std::vector< unsigned int > gradientind; // contains the indices of the gradient images.

		// for each gradient direction ...
		for ( GradientDirectionContainerType::ConstIterator gdcit = this->m_GradientDirectionContainer->Begin(); gdcit != this->m_GradientDirectionContainer->End(); ++gdcit )
		{
			if ( gdcit.Value().one_norm() <= 0.0 )
			{
				baselineind.push_back( gdcit.Index() );
			}
			else
			{
				gradientind.push_back( gdcit.Index() );
			}
		}

		// *************************************************************
		// NEW
		// *************************************************************
		std::cout << "Algorithm: " << m_Algorithm << std::endl;
		std::cout << "Single multi-component image..." << std::endl;

		vnl_vector< double > B( m_NumberOfGradientDirections );
		vnl_vector< double > D( 6 );

		// for each voxel ...
		while ( !git.IsAtEnd() )
		{
			GradientVectorType S = git.Get();

			typename NumericTraits< ReferencePixelType >::AccumulateType S0 = NumericTraits< ReferencePixelType >::Zero;

			// Average the baseline image pixels

			for( unsigned int i = 0; i < baselineind.size(); ++i )
				S0 += S[ baselineind[ i ] ];

			S0 /= this->m_NumberOfBaselineImages;

			TensorPixelType tensor( 0.0 );

			if ( ( S0 != 0 ) && ( S0 >= m_Threshold ) )
			{
				for ( unsigned int i = 0; i < m_NumberOfGradientDirections; i++ )
				{
					if ( S[ gradientind[ i ] ] == 0 )
						B[ i ] = 0;
					else
						B[ i ] = -1 * vcl_log( static_cast< double > ( S[ gradientind[ i ] ] ) / static_cast< double > ( S0 ) ) / this->m_BValue;
				}

				// set in function
				vnl_svd< double > pseudoInverseSolver( m_TensorBasis );

				if ( m_NumberOfGradientDirections > 6 )
					D = pseudoInverseSolver.solve( m_BMatrix * B );
				else
					D = pseudoInverseSolver.solve( B );
				// set in function

				tensor( 0, 0 ) = D[ 0 ];
				tensor( 0, 1 ) = D[ 1 ];
				tensor( 0, 2 ) = D[ 2 ];
				tensor( 1, 1 ) = D[ 3 ];
				tensor( 1, 2 ) = D[ 4 ];
				tensor( 2, 2 ) = D[ 5 ];

				// DEBUG
				std::cout << "Signal: " << S << std::endl;
				std::cout << "B: " << B << std::endl;
				std::cout << "m_BMatrix: " << m_BMatrix << std::endl;
				// DEBUG
			}

			oit.Set( tensor );
			++oit; // Output (reconstructed tensor image) iterator
			++git; // Gradient  image iterator
		}
	}


	/**
	 * Set design matrix.
	 */
	template< class TReferenceImagePixelType, class TGradientImagePixelType, class TTensorPixelType >
	void DiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType, TGradientImagePixelType, TTensorPixelType >::SetDesignMatrix(
			const DesignMatrixType& designMatrix )
	{
		this->m_DesignMatrixContainer = designMatrix;
	}

	/**
	 * Set gradient image.
	 */
	template< class TReferenceImagePixelType, class TGradientImagePixelType, class TTensorPixelType >
	void DiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType, TGradientImagePixelType, TTensorPixelType >::SetGradientImage(
			const GradientImagesType *gradientImage )
	{
		// Make sure crazy users did not call both AddGradientImage and SetGradientImage ...
		if ( m_GradientImageTypeEnumeration == GradientIsInManyImages )
		{
			itkExceptionMacro( << "Cannot call both methods:"
					<< "AddGradientImage and SetGradientImage. Please call only one of them.");
		}

		this->m_NumberOfBaselineImages = 0;

		for ( GradientDirectionContainerType::Iterator it = this->m_GradientDirectionContainer->Begin();
				it != this->m_GradientDirectionContainer->End(); it++ )
		{
			if ( it.Value().one_norm() <= 0.0 )
			{
				this->m_NumberOfBaselineImages++;
			}
			else // Normalize non-zero gradient directions
			{
				it.Value() = it.Value() / it.Value().two_norm();
			}
		}

		this->ProcessObject::SetNthInput( 0, const_cast< GradientImagesType* > ( gradientImage ) );
	}

	/**
	 * Print self.
	 */
	template< class TReferenceImagePixelType, class TGradientImagePixelType, class TTensorPixelType >
	void DiffusionTensor3DReconstructionImageFilter< TReferenceImagePixelType, TGradientImagePixelType, TTensorPixelType >
	::PrintSelf( std::ostream& os, Indent indent ) const
	{
		Superclass::PrintSelf( os, indent );
		os << indent << "Algorithm: " << m_Algorithm << std::endl;
		// TODO
	}

}

#endif
