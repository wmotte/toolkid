#include "tkdCmdParser.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_diag_matrix.h"
#include "vnl/vnl_vector.h"

#include "vnl/algo/vnl_svd.h"

#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "graphCommon.h"

#include "itkImageIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCropImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkLogImageFilter.h"

/**
 * SVD deconvolution (with optional circular deconvolution)
 * to esimate relative perfusion maps (CBF, CBV and MTT).
 *
 * Ostergaard et al, High Resolution Measurement of Cerebral Blood Flow using
 * Intervascular Tracer Bolus Passages. Part I: Mathematical Approach and Statistical Analysis.
 *
 * MRM 36:715-725 (1996)
 *
 * wim@invivonmr.uu.nl, 02-03-2010.
 *
 */
namespace bolus
{
	typedef double PixelType;
	typedef vnl_matrix< PixelType > MatrixType;
	typedef vnl_diag_matrix< PixelType > DiagMatrixType;
	typedef vnl_vector< PixelType > VectorType;

	typedef itk::Image< PixelType, 3 > Image3DType;
	typedef itk::Image< PixelType, 4 > Image4DType;

	typedef itk::ImageLinearIteratorWithIndex< Image4DType > Iterator4DType;
	typedef itk::ImageLinearConstIteratorWithIndex< Image4DType > ConstIterator4DType;

	typedef itk::ImageRegionIteratorWithIndex< Image3DType > Iterator3DType;
	typedef itk::ImageRegionConstIteratorWithIndex< Image3DType > ConstIterator3DType;


	// *******************************************************************************

	/**
	 * Voxel bolus option list.
	 */
	struct voxelBolusParameters
	{
		VectorType cerebralTracer;
		VectorType arterialInputFunction;
		PixelType repetitionTime;
		PixelType toleranceSVD;
		bool circular;
		bool disableFilterAif;
		bool integrateResidueFunction;
	};

	/**
	 * Container for mapping values.
	 */
	struct mapValues
	{
		PixelType cbf;
		PixelType mtt;
		PixelType cbv;
		PixelType chiSquareFit;
		PixelType svdThreshold;
		PixelType oscillationIndex;
		VectorType residue;
	};

	// *******************************************************************************

	/**
	 * Calculate perfusion parameters for voxel time-series.
	 */
	class VoxelBolusSVD
	{

	private:

		/* output parameters */
		mapValues m_MapValues;

	public:

		/**
		 * Run estimation at construction time.
		 */
		VoxelBolusSVD( const voxelBolusParameters& list )
		{
			VectorType aif = list.arterialInputFunction;

			if( list.disableFilterAif )
				aif = PreFilterAif( aif );

			MatrixType A = ConstructA( aif, list.repetitionTime, list.circular );

			MatrixType V;
			DiagMatrixType W;
			MatrixType U;

			Solve( A, V, W, U );

			PixelType svdThreshold = ThresholdW( W, list.toleranceSVD );

			VectorType b = Residuals( V, W, U, list.cerebralTracer );

			// CBV ...
			m_MapValues.cbf = b.max_value();

			// CBF ...
			m_MapValues.cbv = 0;
			for( unsigned int i = 0; i < b.size(); i++ )
			{
				if( list.integrateResidueFunction )
				{
					// integral F_t * R(t) ...
					m_MapValues.cbv += list.repetitionTime * b( i );
				}
				else
				{
					// integral time attenuation curve (rCBV) ..
					m_MapValues.cbv += list.repetitionTime * list.cerebralTracer( i );
				}
			}

			// MTT ...
			m_MapValues.mtt = 0;
			if( m_MapValues.cbf != 0 )
				m_MapValues.mtt = m_MapValues.cbv / m_MapValues.cbf;

			m_MapValues.svdThreshold = svdThreshold;
			m_MapValues.chiSquareFit = -1.0; // TODO
			m_MapValues.oscillationIndex = -1.0; // TODO
			m_MapValues.residue = b;
		}

    /**
     * Threshold accoring Ostergaard.
     */
    PixelType ThresholdW( DiagMatrixType& W, PixelType tolerance )
    {
      VectorType diag = W.diagonal();

      PixelType threshold = diag.max_value() * tolerance;

      for( unsigned int i = 0; i < diag.size(); i++ )
        if( diag( i ) < threshold )
          diag( i ) = 0;

			// insert back ...
			for( unsigned int i = 0; i < diag.size(); i++ )
				W( i, i ) = diag( i );

      return threshold;
    }

		/**
		 * Threshold W, according to percentage tolerance. (Ona Wu).
		 */
		PixelType ThresholdWOld( DiagMatrixType& W, PixelType tolerance )
		{
			VectorType diag = W.diagonal();

			// get S ...
			VectorType S( diag.size(), 0 );
			for( unsigned int i = 0; i < diag.size(); i++ )
				if( diag( i ) != 0 )
					S( i ) = 1 / diag( i );

			// calculate threshold ...
			PixelType threshold = S.max_value() * tolerance;

			// threshold ...
			for( unsigned int i = 0; i < S.size(); i++ )
				if( S( i ) < threshold )
					diag( i ) = 0;

			// insert back ...
			for( unsigned int i = 0; i < diag.size(); i++ )
				W( i, i ) = diag( i );

			return ( 1 / threshold );
		}

		/**
		 * Return output map values.
		 */
		mapValues GetOutput()
		{
			return m_MapValues;
		}

	private:

		/**
		 * Prefilter arterial input function to reduce noise contributions
		 * and to compensate for discretization errors (Check http://boa/wiki/Bolus_tracking).
		 */
		VectorType PreFilterAif( const VectorType& aif )
		{
			VectorType filAif( aif );

			for( unsigned int i = 1; i < filAif.size() - 1; i++ )
				filAif( i ) = ( aif( i - 1 ) + 4 * aif( i ) + aif( i + 1 ) ) / static_cast< PixelType >( 6 );

			return filAif;
		}

		/**
		 * Return deconvolved AIF matrix (Check http://boa/wiki/Bolus_tracking).
		 */
		MatrixType ConstructA( const VectorType& arterialInput, PixelType deltaT, bool circular )
		{
			MatrixType result( arterialInput.size(), arterialInput.size(), 0.0 );

			if( circular ) // TODO -> add zero-padding !
			{
				for( unsigned int r = 0; r < result.rows(); r++ )
					for( unsigned int c = 0; c < result.cols(); c++ )
					{
						if( c >= r )
							result( c, r ) = arterialInput( c - r );
						else
							result( c, r ) = arterialInput( arterialInput.size() - r + c );
					}
			}
			else
			{
				for( unsigned int r = 0; r < result.cols(); r++ )
					for( unsigned int c = r; c < result.rows(); c++ )
					{
						result( c, r ) = arterialInput( c - r );
					}
			}

			return deltaT * result;
		}

		/**
		 * Solve matrix A and return V, W, and U.
		 */
		void Solve( const MatrixType& A, MatrixType& V, DiagMatrixType& W, MatrixType& U )
		{
			vnl_svd< double > svd( A );
			V = svd.V();
			W = svd.W();
			U = svd.U();
		}

		/**
		 * Return residual vector.
		 */
		VectorType Residuals( const MatrixType& V, DiagMatrixType& W, const MatrixType& U, const VectorType& c )
		{
			return V * W * ( U * c );
		}
	};

	// *******************************************************************************

	/**
	 * Bolus option list.
	 */
	struct bolusParameters
	{
		std::string inputFileName;
		std::string maskFileName;
		std::string arterialInputFunctionFileName;
		std::string outputCbvFileName;
		std::string outputCbfFileName;
		std::string outputMttFileName;
		std::string outputRCbvFileName;
		std::string outputChiSquareFileName;
		std::string outputSVDthresholdFileName;
		std::string outputOscillationIndexFileName;

		std::vector< int > range;
		PixelType toleranceSVD;
		PixelType repetitionTime;
		PixelType echoTime;
		bool circular;
		bool disableFilterAif;
		bool integrateResidueFunction;
	};

	// *******************************************************************************

	/**
	 * Calculate perfusion parameters for all voxels within given mask.
	 */
	class BolusSVD
	{

	private:

		Image4DType::Pointer m_rawData;
		Image4DType::Pointer m_RCBVData;
		Image3DType::Pointer m_Mask;
		VectorType m_AIF;

		Image3DType::Pointer m_Cbv;
		Image3DType::Pointer m_Cbf;
		Image3DType::Pointer m_Mtt;
		Image3DType::Pointer m_Chi;
		Image3DType::Pointer m_Oscillation;
		Image3DType::Pointer m_SVDThreshold;

	public:

		/**
		 * Estimate perfusion.
		 */
		void EstimatePerfusion( const bolusParameters& list )
		{
			// [ 1 ] Check input ...
			Checks( list );

			// [ 2 ] Read AIF from file ...
			VectorType AIF;
			graph::Graph< PixelType >::GetVectorFromFile( list.arterialInputFunctionFileName, AIF );

			// [ 3 ] Read data ...
			ReadRaw4D( list.inputFileName );
			ReadMask( list.maskFileName );
			AllocateOutput();

			// [ 4 ] extract relative CBV data and AIF ...
			RCBV( list.echoTime, list.range, AIF );

			// [ 5 ] estimating quantities ...
			unsigned int relativeBolusStartPoint = list.range.at( 2 ) - list.range.at( 0 );
			unsigned int relativeBolusEndPoint = list.range.at( 3 ) - list.range.at( 0 );

			ProcessVoxels( list, relativeBolusStartPoint, relativeBolusEndPoint );

			// [ 6 ] Output maps ...
			Output( list );
		}

	private:

		/**
		 * Output maps, if file names are specified.
		 */
		void Output( const bolusParameters& list )
		{
			if( ! list.outputRCbvFileName.empty() )
			{
				typedef itk::ImageFileWriter< Image4DType > WriterType;
				WriterType::Pointer writer = WriterType::New();
				writer->SetFileName( list.outputRCbvFileName );
				writer->SetInput( m_RCBVData );
				writer->Update();
			}

			if( ! list.outputCbvFileName.empty() )
				Output( list.outputCbvFileName, m_Cbv );

			if( ! list.outputCbfFileName.empty() )
				Output( list.outputCbfFileName, m_Cbf );

			if( ! list.outputMttFileName.empty() )
				Output( list.outputMttFileName, m_Mtt );

			if( ! list.outputChiSquareFileName.empty() )
				Output( list.outputChiSquareFileName, m_Chi );

			if( ! list.outputOscillationIndexFileName.empty() )
				Output( list.outputOscillationIndexFileName, m_Oscillation );

			if( ! list.outputSVDthresholdFileName.empty() )
				Output( list.outputSVDthresholdFileName, m_SVDThreshold );
		}

		/**
		 * Write output.
		 */
		void Output( const std::string& outputFileName, const Image3DType::Pointer& image )
		{
			typedef itk::ImageFileWriter< Image3DType > WriterType;
			WriterType::Pointer writer = WriterType::New();
			writer->SetInput( image );
			writer->SetFileName( outputFileName );

			try
			{
				writer->Update();
			}
			catch( itk::ExceptionObject& e )
			{
				std::cerr << "*** ERROR ***: Could not write output to: " << outputFileName << "!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		/**
		 * Create rCBV data.
		 */
		void RCBV( PixelType echoTime, const std::vector< int >& range, const VectorType& AIF )
		{
			if( echoTime == 0 )
			{
				std::cerr << "*** ERROR ***: echo time should be > 0!" << std::endl;
				exit( EXIT_FAILURE );
			}

			PixelType k = static_cast< PixelType >( 1.0 ) / echoTime;

			unsigned int start = range.at( 0 );
			unsigned int endPre = range.at( 1 );
			unsigned int end = range.at( 3 );
			unsigned int volumes = end - start;

			// crop AIF and convert to rCBV signal...
			m_AIF = VectorType( volumes, 0 	);
			for( unsigned int i = start; i < end; i++ )
				m_AIF( i - start ) = AIF( i );

			PixelType aifPre = 0;
			for( unsigned int i = start; i < endPre; i++ )
				aifPre += AIF( i );

			if( aifPre != 0 )
				aifPre /= ( endPre - start );

			for( unsigned int i = 0; i < m_AIF.size(); i++ )
			{
				if( m_AIF( i ) != 0 )
					m_AIF( i ) = -k * log( m_AIF( i ) / aifPre );
			}

			// extract S(0) volumes ...
			Image4DType::RegionType extractRegion = m_rawData->GetLargestPossibleRegion();
			Image4DType::SizeType extractSize = extractRegion.GetSize();
			Image4DType::IndexType extractIndex = extractRegion.GetIndex();
			extractIndex[3] = start;
			extractSize[3] = endPre - start;
			extractRegion.SetIndex( extractIndex );
			extractRegion.SetSize( extractSize );

			typedef itk::ExtractImageFilter< Image4DType, Image4DType > ExtractImageFilterType;
			ExtractImageFilterType::Pointer cropFilter = ExtractImageFilterType::New();
			cropFilter->SetInput( m_rawData );
			cropFilter->SetExtractionRegion( extractRegion );
			cropFilter->Update();
			Image4DType::Pointer preImages = cropFilter->GetOutput();

			// Calculate average S0 image ...

			Image3DType::Pointer averagePre = Image3DType::New();
			Image3DType::RegionType region;
			Image3DType::SizeType size3D;
			size3D[0] = extractSize[0];
			size3D[1] = extractSize[1];
			size3D[2] = extractSize[2];
			region.SetSize( size3D );
			averagePre->SetRegions( region );
			averagePre->Allocate();
			averagePre->FillBuffer( 0 );

			ConstIterator4DType it( preImages, preImages->GetLargestPossibleRegion() );
			ConstIterator3DType mit( m_Mask, m_Mask->GetLargestPossibleRegion() );
			it.SetDirection( 3 );
			it.GoToBegin();
			mit.GoToBegin();

			unsigned int totalVolumes = extractSize[3];

			while( ! it.IsAtEnd(), ! mit.IsAtEnd() )
			{
				if( mit.Get() != 0 )
				{
					PixelType sum = 0;
					while( ! it.IsAtEndOfLine() )
					{
						sum += it.Get();
						++it;
					}

					Image4DType::IndexType index4D = it.GetIndex();
					Image3DType::IndexType index3D;
					index3D[0] = index4D[0];
					index3D[1] = index4D[1];
					index3D[2] = index4D[2];

					sum /= totalVolumes;
					averagePre->SetPixel( index3D, sum );
				}
				it.NextLine();
				++mit;
			}

			// Crop raw data

			extractSize[3] = end - start;
			extractIndex[3] = start;
			extractRegion.SetSize( extractSize );
			extractRegion.SetIndex( extractIndex );
			ExtractImageFilterType::Pointer newCrop = ExtractImageFilterType::New();
			newCrop->SetExtractionRegion( extractRegion );
			newCrop->SetInput( m_rawData );
			newCrop->Update();
			m_rawData = newCrop->GetOutput();

			// Convert data to rCBV ...

			Image4DType::SizeType alloSize = m_rawData->GetLargestPossibleRegion().GetSize();
			Image4DType::IndexType alloIndex = m_rawData->GetLargestPossibleRegion().GetIndex();
			Image4DType::RegionType alloRegion;
			alloRegion.SetSize( alloSize );
			alloRegion.SetIndex( alloIndex );
			m_RCBVData = Image4DType::New();
			m_RCBVData->SetRegions( alloRegion );
			m_RCBVData->Allocate();
			m_RCBVData->FillBuffer( 0 );

			ConstIterator4DType iit( m_rawData, m_rawData->GetLargestPossibleRegion() );
			Iterator4DType rit( m_RCBVData, m_RCBVData->GetLargestPossibleRegion() );
			ConstIterator3DType miit( m_Mask, m_Mask->GetLargestPossibleRegion() );

			iit.SetDirection( 3 );
			rit.SetDirection( 3 );
			iit.GoToBegin();
			rit.GoToBegin();
			miit.GoToBegin();

			while( ! iit.IsAtEnd() )
			{
				if( miit.Get() != 0 )
				{
					Image4DType::IndexType index4D = iit.GetIndex();
					Image3DType::IndexType index3D;
					index3D[0] = index4D[0];
					index3D[1] = index4D[1];
					index3D[2] = index4D[2];

					PixelType preValue = averagePre->GetPixel( index3D );

					while( ! iit.IsAtEndOfLine() )
					{
						if( preValue != 0 )
						{
							rit.Set( -k * ( log( iit.Get() / preValue ) ) );
							++iit;
							++rit;
						}
					}
				}
				iit.NextLine();
				rit.NextLine();
				++miit;
			}
		}

		/**
		 * Read 3D volume.
		 */
		void ReadMask( const std::string& inputFileName )
		{
			if( ! inputFileName.empty() )
			{
				typedef itk::ImageFileReader< Image3DType > ReaderType;
				ReaderType::Pointer reader = ReaderType::New();
				reader->SetFileName( inputFileName );
				reader->Update();
				m_Mask = reader->GetOutput();
			}
			else
			{
				// use 4D input data as reference for full mask...
				Image4DType::SizeType size4D = m_rawData->GetLargestPossibleRegion().GetSize();
				Image4DType::IndexType index4D = m_rawData->GetLargestPossibleRegion().GetIndex();
				Image3DType::SizeType size3D;
				Image3DType::IndexType index3D;
				size3D[0] = size4D[0];
				size3D[1] = size4D[1];
				size3D[2] = size4D[2];
				index3D[0] = index4D[0];
				index3D[1] = index4D[1];
				index3D[2] = index4D[2];
				Image3DType::RegionType region3D;
				region3D.SetSize( size3D );
				region3D.SetIndex( index3D );
				m_Mask = Image3DType::New();
				m_Mask->SetRegions( region3D );
				m_Mask->Allocate();
				m_Mask->FillBuffer( 1 );
			}
		}

		/**
		 * Allocate space for output 3D maps.
		 */
		void AllocateOutput()
		{
			// use 4D input data as reference for full mask...
			Image4DType::SizeType size4D = m_rawData->GetLargestPossibleRegion().GetSize();
			Image4DType::IndexType index4D = m_rawData->GetLargestPossibleRegion().GetIndex();
			Image3DType::SizeType size3D;
			Image3DType::IndexType index3D;
			size3D[0] = size4D[0];
			size3D[1] = size4D[1];
			size3D[2] = size4D[2];
			index3D[0] = index4D[0];
			index3D[1] = index4D[1];
			index3D[2] = index4D[2];
			Image3DType::RegionType region3D;
			region3D.SetSize( size3D );
			region3D.SetIndex( index3D );

			Image4DType::SpacingType spacing4D = m_rawData->GetSpacing();
			Image3DType::SpacingType spacing3D;
			spacing3D[0] = spacing4D[0];
			spacing3D[1] = spacing4D[1];
			spacing3D[2] = spacing4D[2];

			Image4DType::PointType origin4D = m_rawData->GetOrigin();
			Image3DType::PointType origin3D;
			origin3D[0] = origin4D[0];
			origin3D[1] = origin4D[1];
			origin3D[2] = origin4D[2];

			m_Cbv = Image3DType::New();
			m_Cbv->SetRegions( region3D );
			m_Cbv->SetSpacing( spacing3D );
			m_Cbv->SetOrigin( origin3D );
			m_Cbv->Allocate();
			m_Cbv->FillBuffer( 0 );

			m_Cbf = Image3DType::New();
			m_Cbf->SetRegions( region3D );
			m_Cbf->SetSpacing( spacing3D );
			m_Cbf->SetOrigin( origin3D );
			m_Cbf->Allocate();
			m_Cbf->FillBuffer( 0 );

			m_Mtt = Image3DType::New();
			m_Mtt->SetRegions( region3D );
			m_Mtt->SetSpacing( spacing3D );
			m_Mtt->SetOrigin( origin3D );
			m_Mtt->Allocate();
			m_Mtt->FillBuffer( 0 );

			m_Chi = Image3DType::New();
			m_Chi->SetRegions( region3D );
			m_Chi->SetSpacing( spacing3D );
			m_Chi->SetOrigin( origin3D );
			m_Chi->Allocate();
			m_Chi->FillBuffer( 0 );

			m_Oscillation = Image3DType::New();
			m_Oscillation->SetRegions( region3D );
			m_Oscillation->SetSpacing( spacing3D );
			m_Oscillation->SetOrigin( origin3D );
			m_Oscillation->Allocate();
			m_Oscillation->FillBuffer( 0 );

			m_SVDThreshold = Image3DType::New();
			m_SVDThreshold->SetRegions( region3D );
			m_SVDThreshold->SetSpacing( spacing3D );
			m_SVDThreshold->SetOrigin( origin3D );
			m_SVDThreshold->Allocate();
			m_SVDThreshold->FillBuffer( 0 );
		}

		/**
		 * Read raw data.
		 */
		void ReadRaw4D( const std::string& inputFileName )
		{
			typedef itk::ImageFileReader< Image4DType > ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName( inputFileName );
			reader->Update();
			m_rawData = reader->GetOutput();
		}

		/**
		 * Process all voxels within mask.
		 */
		void ProcessVoxels( const bolusParameters& list, unsigned int startIndex, unsigned int endIndex )
		{
			// set aif to given indices...

			VectorType aif( endIndex - startIndex, 0 );
			for( unsigned int i = 0; i < aif.size(); i++ )
				aif( i ) = m_AIF( i + startIndex );

			// set rcbv region to given indices...

			Image4DType::RegionType rcbvRegion = m_RCBVData->GetLargestPossibleRegion();
			Image4DType::SizeType rcbvSize = rcbvRegion.GetSize();
			Image4DType::IndexType rcbvIndex = rcbvRegion.GetIndex();
			rcbvSize[3] = endIndex - startIndex;
			rcbvIndex[3] = startIndex;
			rcbvRegion.SetSize( rcbvSize );
			rcbvRegion.SetIndex( rcbvIndex );

			ConstIterator4DType it( m_RCBVData, rcbvRegion );
			ConstIterator3DType mit( m_Mask, m_Mask->GetLargestPossibleRegion() );

			it.SetDirection( 3 );
			it.GoToBegin();
			mit.GoToBegin();

			while( ! it.IsAtEnd(), ! mit.IsAtEnd() )
			{
				if( mit.Get() != 0 )
				{
					VectorType cerebralTracer( rcbvSize[3], 0 );

					unsigned int i = 0;
					while( ! it.IsAtEndOfLine() )
					{
						cerebralTracer( i ) = it.Get();
						++it;
						++i;
					}

					voxelBolusParameters p;
					p.cerebralTracer = cerebralTracer;
					p.arterialInputFunction = aif;
					p.repetitionTime = list.repetitionTime;
					p.circular = list.circular;
					p.disableFilterAif = list.disableFilterAif;
					p.toleranceSVD = list.toleranceSVD;
					p.integrateResidueFunction = list.integrateResidueFunction;
					VoxelBolusSVD svd( p );

					Image4DType::IndexType index4D = it.GetIndex();
					Image3DType::IndexType index3D;
					index3D[0] = index4D[0];
					index3D[1] = index4D[1];
					index3D[2] = index4D[2];

					FillVoxel( svd.GetOutput(), index3D );
				}
				it.NextLine();
				++mit;
			}
		}

		/**
		 * Fill output files with values.
		 */
		void FillVoxel( const mapValues& mapValue, const Image3DType::IndexType& index )
		{
			m_Cbv->SetPixel( index, mapValue.cbv );
			m_Cbf->SetPixel( index, mapValue.cbf );
			m_Mtt->SetPixel( index, mapValue.mtt );
			m_Chi->SetPixel( index, mapValue.chiSquareFit );
			m_Oscillation->SetPixel( index, mapValue.oscillationIndex );
			m_SVDThreshold->SetPixel( index, mapValue.svdThreshold );
		}

		/**
		 * Check parameters, and exist if not correct.
		 */
		void Checks( const bolusParameters& list )
		{
			// input ...
			if ( graph::Graph< PixelType >::GetImageDimensions( list.inputFileName ) != 4 )
			{
				std::cout << "*** ERROR ***: Input \"" << list.inputFileName << "\" is not a 4D image!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// mask ...
			if( ! list.maskFileName.empty() )
			{
				if ( graph::Graph< PixelType >::GetImageDimensions( list.maskFileName ) != 3 )
				{
					std::cout << "*** ERROR ***: Mask \"" << list.maskFileName << "\" is not a 3D image!" << std::endl;
					exit( EXIT_FAILURE );
				}
			}

			// tolerance ...
			if( list.toleranceSVD < 0 || list.toleranceSVD > 1 )
			{
				std::cout << "*** ERROR ***: SVD tolerance should have a value between 0 and 1!" << std::endl;
				exit( EXIT_FAILURE );
			}

			if( ( list.range.empty() ) || ( list.range.size() != 4 ) )
			{
				std::cerr << "*** ERROR ***: Range should contain 4 values! (start pre, end pre, bolus, end volume)" << std::endl;
				exit( EXIT_FAILURE );
			}

			if( ! ( list.range.at( 0 ) < list.range.at( 1 ) ) )
			{
				std::cerr << "*** ERROR ***: start pre exceeds end pre volume!" << std::endl;
				std::cerr << "[Valid range : 0 < start pre < end pre < bolus < end volume]" << std::endl;
				exit( EXIT_FAILURE );
			}

			if( ! ( list.range.at( 1 ) < list.range.at( 2 ) ) )
			{
				std::cerr << "*** ERROR ***: end pre volume exceeds bolus volume!" << std::endl;
				std::cerr << "[Valid range : 0 < start pre < end pre < bolus < end volume]" << std::endl;
				exit( EXIT_FAILURE );
			}

			if( ! ( list.range.at( 2 ) < list.range.at( 3 ) ) )
				{
					std::cerr << "*** ERROR ***: bolus volume exceeds end volume!" << std::endl;
					std::cerr << "[Valid range : 0 < start volume < bolus volume < end volume]" << std::endl;
					exit( EXIT_FAILURE );
				}
		}
	};



} // end namespace bolus

/**
 * Main.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Singular value decomposition method for dynamic susceptibility contrast-weighted (DSC) MR images." );

	bolus::bolusParameters args;

	args.toleranceSVD = 0.2;
	args.repetitionTime = 0.5;
	args.echoTime = 0.019;
	args.circular = false;
	args.disableFilterAif = false;
	args.integrateResidueFunction = false;


	p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" ) ->SetInput( "string" ) ->SetDescription( "Perfusion 4D input image" )->SetMinMax( 1,1 )->SetRequired( true );
	p.AddArgument( args.arterialInputFunctionFileName, "arterial-input-function" ) ->AddAlias( "a" ) ->SetInput( "string" ) ->SetDescription( "Arterial input function (column format)" )->SetMinMax( 1,1 )->SetRequired( true );
	p.AddArgument( args.range, "indices" ) ->AddAlias( "r" ) ->SetInput( "ints" ) ->SetDescription( "(1) Start pre, (2) End pre, (3) Bolus, (4) End volume" )->SetMinMax( 4, 4 )->SetRequired( true );
	p.AddArgument( args.maskFileName, "mask" ) ->AddAlias( "m" ) ->SetInput( "string" ) ->SetDescription( "Mask 3D image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputCbvFileName, "output-cbv" ) ->AddAlias( "cbv" ) ->SetInput( "string" ) ->SetDescription( "Output CBV image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputCbfFileName, "output-cbf" ) ->AddAlias( "cbf" ) ->SetInput( "string" ) ->SetDescription( "Output CBF image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputMttFileName, "output-mtt" ) ->AddAlias( "mtt" ) ->SetInput( "string" ) ->SetDescription( "Output MTT image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputRCbvFileName, "output-rcbv" ) ->AddAlias( "rcbv" ) ->SetInput( "string" ) ->SetDescription( "Output 4D rCBV image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputChiSquareFileName, "output-chi-square" ) ->AddAlias( "chi" ) ->SetInput( "string" ) ->SetDescription( "Output Chi square goodness of fit image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputSVDthresholdFileName, "output-svd-threshold" ) ->AddAlias( "svdt" ) ->SetInput( "string" ) ->SetDescription( "Output SVD threshold image" )->SetMinMax( 1,1 );
	p.AddArgument( args.outputOscillationIndexFileName, "output-oscillation-index" ) ->AddAlias( "oi" ) ->SetInput( "string" ) ->SetDescription( "Output oscillation index image" )->SetMinMax( 1,1 );
	p.AddArgument( args.toleranceSVD, "tolerance-svd" ) ->AddAlias( "tsvd" ) ->SetInput( "float" ) ->SetDescription( "SVD elimination threshold diagonal elements W matrix (default: 0.2 ~ 20%)" )->SetMinMax( 1,1 );
	p.AddArgument( args.repetitionTime, "repetition-time" ) ->AddAlias( "tr" ) ->SetInput( "float" ) ->SetDescription( "Repetition time in seconds (default: 0.5)" )->SetMinMax( 1, 1 );
	p.AddArgument( args.echoTime, "echo-time" ) ->AddAlias( "te" ) ->SetInput( "float" ) ->SetDescription( "Echo time in seconds (default: 0.019)" )->SetMinMax( 1, 1 );
	p.AddArgument( args.circular, "circular" ) ->AddAlias( "c" ) ->SetInput( "bool" ) ->SetDescription( "Circular deconvolution (default: false)" );
	p.AddArgument( args.disableFilterAif, "disable-filter-aif" ) ->AddAlias( "f" ) ->SetInput( "bool" ) ->SetDescription( "Disable prefilter arterial input function (default: false)" );
	p.AddArgument( args.integrateResidueFunction, "integrate-residue" ) ->AddAlias( "ir" ) ->SetInput( "bool" ) ->SetDescription( "Integrate F_t * R(t) to estimate CBV (default: integrate time attenuation curve ~ relative CBV)" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	bolus::BolusSVD bolusSvd;

	bolusSvd.EstimatePerfusion( args );

	return EXIT_SUCCESS;
}


