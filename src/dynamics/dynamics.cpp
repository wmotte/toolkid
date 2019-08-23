#include "tkdCmdParser.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkExtractImageFilter.h"

#include <vnl/algo/vnl_lsqr.h>
#include <vnl/vnl_sparse_matrix_linear_system.h>
#include <vnl/vnl_least_squares_function.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_sparse_matrix.h"

/**
 * wim@invivonmr.uu.nl, 30-08-2010.
 *
 * --------------
 * Patlak method:
 * --------------
 * C. S. Patlak, R. G. Blasberg (1985).
 * "Graphical evaluation of blood-to-brain transfer constants from multiple-time uptake data. Generalizations".
 * Journal of Cerebral Blood Flow and Metabolism 5 (4): 584–590.
 *
 * --------------
 * Ichise method:
 * --------------
 * M Ichise, H Toyama, RB Innis, RE Carson (2002)
 * "Strategies to improve neuroreceptor parameter estimation by linear regression analysis".
 * Journal of Cerebral Blood Flow and Metabolism 22 (10): 1271–1281.
 *
 */
namespace dynamics
{
	typedef double PixelType;

    typedef vnl_vector< PixelType > VectorType;
    typedef vnl_sparse_matrix< PixelType > SparseMatrixType;
	
    typedef itk::Image< PixelType, 3 > Image3DType;
	typedef itk::ImageFileReader< Image3DType > Image3DFileReaderType;
	typedef itk::ImageFileWriter< Image3DType > Image3DFileWriterType;
	
	typedef itk::Image< PixelType, 4 > Image4DType;
    typedef itk::ImageFileReader< Image4DType > Image4DFileReaderType;

    typedef itk::ImageLinearIteratorWithIndex< Image4DType > Iterator4DType;
    typedef itk::ImageLinearConstIteratorWithIndex< Image4DType > ConstIterator4DType;

    typedef itk::ImageRegionIteratorWithIndex< Image3DType > Iterator3DType;
    typedef itk::ImageRegionConstIteratorWithIndex< Image3DType > ConstIterator3DType;

    typedef itk::ExtractImageFilter< Image4DType, Image4DType > ExtractImageFilterType;

	// *******************************************************************************

	/**
	 * Option list.
	 */
	struct parameters
	{
		std::string inputFileName;
		std::string maskFileName;
		std::string sinusFileName;
		std::string outputVFileName;
		std::string outputKFileName;
		std::string outputResidualsFileName;
		double hematocrit;
		double relaxivity;
		double TR;
		int method;
		bool verbose;
	};

	// *******************************************************************************

	/**
	 * 1 => Patlak (irreversible compartment, ordinary least squares)
	 * 2 => Ichise (irreversible compartment, multilinear least squares)
	 * 3 => Average of 1 and 2
	 */
	class GraphicalAnalysis
	{

	public:

		/**
		 * Run Patlak.
		 */
		void Run( parameters& args )
		{
			Checks( args );

			if( args.verbose )
				std::cout << "Create concentration maps..." << std::endl;

			Image4DType::Pointer deltaR1 = GetConcentrationMaps( args.inputFileName, args.relaxivity );

			if( args.verbose )
				std::cout << "Read sinus and/or mask image..." << std::endl;

			Image3DType::Pointer sinus = ReadImage3D( args.sinusFileName );
		    Image3DType::Pointer mask = GetMask( sinus, args.maskFileName );

		    if( args.verbose )
		    	std::cout << "Get average sinus signal and integral..." << std::endl;

		    VectorType Ct = AverageSinus( deltaR1, sinus, args.hematocrit );
			VectorType Ct_int = Integral( Ct, args.TR );

			Image3DType::Pointer interceptImage = Image3DType::New();
			Image3DType::Pointer slopeImage = Image3DType::New();
			Image3DType::Pointer residualImage = Image3DType::New();

	    	if( args.method == 1 ) // Patlak 1985
	    	{
	    		if( args.verbose )
	    			std::cout << "Fit Patlak model..." << std::endl;

	    		FitPatlak( deltaR1, mask, Ct, Ct_int, interceptImage, slopeImage, residualImage );
	    	}
	    	else if( args.method == 2 ) // Ichise 2002
			{
	    		if( args.verbose )
	    			std::cout << "Fit Ichise model..." << std::endl;

	    		FitIchise( deltaR1, mask, Ct, Ct_int, interceptImage, slopeImage, residualImage );
			}
	    	else if( args.method == 3 ) // Average Patlak & Ichise (see conclusion review Logan 2003)
	    	{
	    		Image3DType::Pointer interceptImage2 = Image3DType::New();
	    		Image3DType::Pointer slopeImage2 = Image3DType::New();
	    		Image3DType::Pointer residualImage2 = Image3DType::New();

	    		if( args.verbose )
	    			std::cout << "Fit Ichise model..." << std::endl;

	    		FitIchise( deltaR1, mask, Ct, Ct_int, interceptImage, slopeImage, residualImage );

	    		if( args.verbose )
	    			std::cout << "Fit Patlak model..." << std::endl;

	    		FitPatlak( deltaR1, mask, Ct, Ct_int, interceptImage2, slopeImage2, residualImage2 );

	    		if( args.verbose )
	    			std::cout << "Average Patlak and Ichise maps..." << std::endl;

	    		interceptImage = Average( interceptImage, interceptImage2 );
	    		slopeImage = Average( slopeImage, slopeImage2 );
	    		residualImage = Average( residualImage, residualImage2 );
	    	}
	    	else
	    	{
	    		std::cout << "*** ERROR ***: method not implemented (yet)!" << std::endl;
	    		exit( EXIT_FAILURE );
	    	}

	    	if( args.verbose )
	    		std::cout << "Write output..." << std::endl;

	    	// Write
	    	if( ! args.outputKFileName.empty() )
				WriteImage( interceptImage, args.outputKFileName );
			if( ! args.outputVFileName.empty() )
				WriteImage( slopeImage, args.outputVFileName );
			if( ! args.outputResidualsFileName.empty() )
				WriteImage( residualImage, args.outputResidualsFileName );

	    	exit( EXIT_SUCCESS );
		}

	private:

		/**
		 * Average images.
		 */
		Image3DType::Pointer Average( const Image3DType::Pointer& V1, const Image3DType::Pointer& V2 )
		{

			Image3DType::Pointer output = Image3DType::New();
			output->CopyInformation( V1 );
			output->SetRegions( V1->GetLargestPossibleRegion() );
			output->Allocate();
			output->FillBuffer( 0 );

			ConstIterator3DType it1( V1, V1->GetLargestPossibleRegion() );
			ConstIterator3DType it2( V2, V2->GetLargestPossibleRegion() );
			Iterator3DType out( output, output->GetLargestPossibleRegion() );

			it1.GoToBegin();
			it2.GoToBegin();
			out.GoToBegin();

			while( ! it1.IsAtEnd(), ! it2.IsAtEnd(), ! out.IsAtEnd() )
			{
				out.Set( ( it1.Get() + it2.Get() ) / 2 );
				++it1;
				++it2;
				++out;
			}

			return output;
		}

		/**
		 * Ichise fit.
		 *
		 * Check: JCBFM 2002.
		 */
		void FitIchise( const Image4DType::Pointer R1, const Image3DType::Pointer& mask,
				const VectorType Ct, const VectorType Cpt_int,
				Image3DType::Pointer& interceptImage, Image3DType::Pointer& slopeImage, Image3DType::Pointer& residualImage )
		{
			Image4DType::RegionType region4D = R1->GetLargestPossibleRegion();
			Image4DType::PointType origin4D = R1->GetOrigin();
			Image4DType::SpacingType spacing4D = R1->GetSpacing();

			unsigned int totalVolumes = region4D.GetSize()[3];

			if( Ct.size() != totalVolumes )
			{
				std::cerr << "*** ERROR ***: Ct size different from R1 size!" << std::endl;
				exit( EXIT_FAILURE );
			}
			if( Ct.size() != Cpt_int.size() )
			{
				std::cerr << "*** ERROR ***: Ct size different from Cpt_int size!" << std::endl;
				exit( EXIT_FAILURE );
			}

			Image3DType::RegionType region = mask->GetLargestPossibleRegion();
			Image3DType::PointType origin;
			Image3DType::SpacingType spacing;

			origin[0] = origin4D[0];
			origin[1] = origin4D[1];
			origin[2] = origin4D[2];

			spacing[0] = spacing4D[0];
			spacing[1] = spacing4D[1];
			spacing[2] = spacing4D[2];

			interceptImage->SetRegions( region );
			interceptImage->Allocate();
			interceptImage->FillBuffer( 0 );
			interceptImage->SetOrigin( origin );
			interceptImage->SetSpacing( spacing );

			slopeImage->SetRegions( region );
			slopeImage->Allocate();
			slopeImage->FillBuffer( 0 );
			slopeImage->SetOrigin( origin );
			slopeImage->SetSpacing( spacing );

			residualImage->SetRegions( region );
			residualImage->Allocate();
			residualImage->FillBuffer( 0 );
			residualImage->SetOrigin( origin );
			residualImage->SetSpacing( spacing );

			ConstIterator4DType it( R1, R1->GetLargestPossibleRegion() );
			ConstIterator3DType mit( mask, mask->GetLargestPossibleRegion() );
			it.SetDirection( 3 );
			it.GoToBegin();
			mit.GoToBegin();

			while( ! it.IsAtEnd(), ! mit.IsAtEnd() )
			{
				if( mit.Get() != 0 )
				{
					// init input matrix and intercept
					SparseMatrixType X( totalVolumes, 2 );

					VectorType V1( totalVolumes );
					VectorType V2( totalVolumes );
					V1.fill( 0 );
					V2.fill( 0 );

					while( ! it.IsAtEndOfLine() )
					{
						V1( it.GetIndex()[3] ) = Cpt_int( it.GetIndex()[3] );
						V2( it.GetIndex()[3] ) = it.Get();
						++it;
					}

					for( unsigned int i = 0; i < V1.size(); i++ )
					{
						X( i, 0 ) = V1( i );
						X( i, 1 ) = V2( i );
					}

					// Fit lm( dependent ~ X )
					VectorType coef = FitLine( Ct, X );

					Image4DType::IndexType index4D = it.GetIndex();
					Image3DType::IndexType index3D;
					index3D[0] = index4D[0];
					index3D[1] = index4D[1];
					index3D[2] = index4D[2];

					interceptImage->SetPixel( index3D, 1.0 / coef( 1 ) ); // K (1/b)
					slopeImage->SetPixel( index3D, -1 * ( coef( 0 ) / coef( 1 ) ) ); // V (-beta1/beta2)
					residualImage->SetPixel( index3D, coef( 2 ) );
				}
				it.NextLine();
				++mit;
			}
		}

		/**
		 * Multilinear fit. TODO, generalize for more than 2 columns.
		 */
		VectorType FitLine( const VectorType& y, const SparseMatrixType& X )
		{
			vnl_sparse_matrix_linear_system< PixelType > ls( X, y );
			vnl_lsqr lsqr( ls );
			VectorType b( 2 );
			b[0] = b[1] = 0.0;

			// get coefficients ...
			lsqr.minimize( b );

			//lsqr.diagnose_outcome( std::cerr );

			VectorType output( 3 );
			output( 0 ) = b[ 0 ]; // intercept
			output( 1 ) = b[ 1 ]; // slope
			output( 2 ) = lsqr.get_resid_norm_estimate();

			return output;
		}

		/**
		 * Read mask, if specified, else create image with all voxels set to 1.
		 */
		Image3DType::Pointer GetMask( const Image3DType::Pointer& example, const std::string& maskFileName )
		{
			Image3DType::Pointer mask = Image3DType::New();

			if( maskFileName.size() > 0 )
				mask = ReadImage3D( maskFileName );
			else
			{
				mask->CopyInformation( example );
				mask->SetRegions( example->GetLargestPossibleRegion() );
				mask->Allocate();
				mask->FillBuffer( 1.0 );
			}

			return mask;
		}

		/**
		 * Check input parameters and files.
		 */
		void Checks( const parameters& args )
		{
			if( args.outputKFileName.empty() && args.outputResidualsFileName.empty() && args.outputVFileName.empty() )
			{
				std::cout << "*** ERROR ***: no output specified!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// Check dims.
			if( ( GetDimensions( args.inputFileName ) != 4 )
				|| ( GetDimensions( args.sinusFileName ) != 3 ) )
			{
				std::cerr << "*** ERROR ***: Check input/mask/sinus image dimensions!" << std::endl;
				exit( EXIT_FAILURE );
			}

			if( ( args.hematocrit < 0 ) || ( args.hematocrit > 1 ) )
			{
				std::cerr << "*** ERROR ***: hematocrit out of bounds!" << std::endl;
				exit( EXIT_FAILURE );
			}

			if( args.relaxivity <= 0 )
			{
				std::cerr << "*** ERROR ***: Gd-DTPA relaxivity should be > 0!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		/**
		 * Construct Cint(t) = Integral{ C(t)dt} / C(t)
		 */
		VectorType Integral( const VectorType C, double TR )
		{
			VectorType Cint( C.size() );
			Cint.fill( 0 );

			for( unsigned int t = 0; t < C.size(); t++ )
			{
				double Ci = 0;

				for( unsigned int i = 0; i < t; i++ )
					Ci += C( i ) * TR;

				Cint( t ) = Ci;
			}

			return Cint;
		}

		/**
		 * Get average signal from 4D given mask.
		 */
		VectorType AverageSinus( const Image4DType::Pointer image, const Image3DType::Pointer& mask, double hematocrit )
		{
			ConstIterator4DType it( image, image->GetLargestPossibleRegion() );
			ConstIterator3DType mit( mask, mask->GetLargestPossibleRegion() );
			it.SetDirection( 3 );
			it.GoToBegin();
			mit.GoToBegin();

			Image4DType::RegionType region = image->GetLargestPossibleRegion();

			unsigned int totalVolumes = region.GetSize()[3];

			VectorType output( totalVolumes );
			output.fill( 0 );

			unsigned int nonZeroVoxels = 0;

			while( ! it.IsAtEnd(), ! mit.IsAtEnd() )
			{
				if( mit.Get() != 0 )
				{
					while( ! it.IsAtEndOfLine() )
					{
						output( it.GetIndex()[3] ) += it.Get();
						++it;
					}
					++nonZeroVoxels;
				}
				it.NextLine();
				++mit;
			}

			for( unsigned int i = 0; i < output.size(); i++ )
			{
				output( i ) /= nonZeroVoxels;
				output( i ) /= ( 1 - hematocrit );
			}

			return output;
		}

		/**
		 * Patlak fit.
		 *
		 * Check: http://en.wikipedia.org/wiki/Patlak_plot
		 */
		void FitPatlak( const Image4DType::Pointer R1,
				  const Image3DType::Pointer& mask,
				  const VectorType C,
				  const VectorType Cint,
				Image3DType::Pointer& interceptImage,
				Image3DType::Pointer& slopeImage,
				Image3DType::Pointer& residualImage )
		{
			Image4DType::RegionType region4D = R1->GetLargestPossibleRegion();
			Image4DType::PointType origin4D = R1->GetOrigin();
			Image4DType::SpacingType spacing4D = R1->GetSpacing();

			unsigned int totalVolumes = region4D.GetSize()[3];

            if( C.size() != totalVolumes )
            {
            	std::cerr << "*** ERROR ***: C size different from R1 size!" << std::endl;
            	exit( EXIT_FAILURE );
            }
            if( C.size() != Cint.size() )
			{
				std::cerr << "*** ERROR ***: C size different from Cint size!" << std::endl;
				exit( EXIT_FAILURE );
			}

			Image3DType::RegionType region = mask->GetLargestPossibleRegion();
			Image3DType::PointType origin;
			Image3DType::SpacingType spacing;

			origin[0] = origin4D[0];
			origin[1] = origin4D[1];
			origin[2] = origin4D[2];

			spacing[0] = spacing4D[0];
			spacing[1] = spacing4D[1];
			spacing[2] = spacing4D[2];

            interceptImage->SetRegions( region );
            interceptImage->Allocate();
            interceptImage->FillBuffer( 0 );
            interceptImage->SetOrigin( origin );
            interceptImage->SetSpacing( spacing );

            slopeImage->SetRegions( region );
            slopeImage->Allocate();
            slopeImage->FillBuffer( 0 );
            slopeImage->SetOrigin( origin );
            slopeImage->SetSpacing( spacing );

            residualImage->SetRegions( region );
            residualImage->Allocate();
            residualImage->FillBuffer( 0 );
            residualImage->SetOrigin( origin );
            residualImage->SetSpacing( spacing );

            ConstIterator4DType it( R1, R1->GetLargestPossibleRegion() );
            ConstIterator3DType mit( mask, mask->GetLargestPossibleRegion() );
            it.SetDirection( 3 );
            it.GoToBegin();
            mit.GoToBegin();

            while( ! it.IsAtEnd(), ! mit.IsAtEnd() )
            {
                if( mit.Get() != 0 )
                {
                    VectorType Rt( totalVolumes );
                    Rt.fill( 0 );

                    while( ! it.IsAtEndOfLine() )
                    {
                        Rt( it.GetIndex()[3] ) = it.Get() / C( it.GetIndex()[3] );
                        ++it;
                    }

                    // Fit lm( Rt ~ Cint )
                    VectorType coef = FitLine( Rt, Cint );

                    Image4DType::IndexType index4D = it.GetIndex();
                    Image3DType::IndexType index3D;
                    index3D[0] = index4D[0];
                    index3D[1] = index4D[1];
                    index3D[2] = index4D[2];

                    interceptImage->SetPixel( index3D, coef( 0 ) );
                    slopeImage->SetPixel( index3D, coef( 1 ) );
                    residualImage->SetPixel( index3D, coef( 2 ) );
                }
                it.NextLine();
                ++mit;
            }
		}

		/**
		 * Check dimensions of inputfile. In case of error, EXIT_FAILURE is returned.
		 */
		unsigned int GetDimensions( const std::string& inputFileName )
		{
			itk::ImageIOFactory::ImageIOBasePointer io =
					itk::ImageIOFactory::CreateImageIO( inputFileName.c_str(), itk::ImageIOFactory::ReadMode );
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
		 * Least squares fit => y ~ 1 + x (intercept is always added)
		 *
		 * fit linear regression line. y = b * x + a
		 * return vector with:
		 * [0] => a
		 * [1] => b
		 * [2] => residual norm
		 */
		VectorType FitLine( const VectorType& y, const VectorType& x )
		{

			// check...
			if ( x.size() != y.size() )
			{
				std::cout << "*** ERROR ***: FitLine, inputs do not match in size!" << std::endl;
				exit( EXIT_FAILURE );
			}

			// init input matrix and intercept
			SparseMatrixType X( x.size(), 2 );

			for( unsigned int i = 0; i < x.size(); i++ )
			{
				X( i, 0 ) = 1;
				X( i, 1 ) = x( i );
			}

			vnl_sparse_matrix_linear_system< PixelType > ls( X, y );
			vnl_lsqr lsqr( ls );
			VectorType b( 2 );
			b[0] = b[1] = 0.0;

			// get coefficients ...
			lsqr.minimize( b );

			//lsqr.diagnose_outcome( std::cerr );

			VectorType output( 3 );
			output( 0 ) = b[ 0 ]; // intercept
			output( 1 ) = b[ 1 ]; // slope
			output( 2 ) = lsqr.get_resid_norm_estimate();

			return output;
		}

		/**
		 * Write.
		 */
		void WriteImage( const Image3DType::Pointer image, const std::string outputFileName )
		{
			Image3DFileWriterType::Pointer writer = Image3DFileWriterType::New();
			writer->SetFileName( outputFileName );
			writer->SetInput( image );

			try
			{
				writer->Update();
			} catch( ... )
			{
				std::cout << "*** ERROR ***: unable to write to " << outputFileName << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		/**
		 * Read 3d image.
		 */
		Image3DType::Pointer ReadImage3D( const std::string& image )
		{
			Image3DFileReaderType::Pointer reader = Image3DFileReaderType::New();
			reader->SetFileName( image );
			reader->Update();
			return reader->GetOutput();
		}
		
        /**
		 * Read 4d image and create Delta R1 images.
		 */
		Image4DType::Pointer GetConcentrationMaps( const std::string& imageName, double relaxivity )
		{
			Image4DFileReaderType::Pointer reader = Image4DFileReaderType::New();
			reader->SetFileName( imageName );
			reader->Update();

			// init ratio
			Image4DType::Pointer image = reader->GetOutput();
			Image4DType::Pointer ratio = Image4DType::New();
			ratio->CopyInformation( image );
			ratio->SetRegions( image->GetLargestPossibleRegion() );
			ratio->Allocate();
			ratio->FillBuffer( 0 );

			// iterator
			ConstIterator4DType it( image, image->GetLargestPossibleRegion() );
			it.SetDirection( 3 );
			it.GoToBegin();

			while( ! it.IsAtEnd() )
			{
				// divide all pixels with pixel value at t = 0
				PixelType first = it.Get();
				while( ! it.IsAtEndOfLine() )
				{
					// C(t) = ( R1(t) - R1(0) ) / relaxivity
					ratio->SetPixel( it.GetIndex(), ( it.Get() - first ) / relaxivity );
					++it;
				}
				it.NextLine();
			}

			return RemoveFirstImage( ratio );
		}

		/**
		 * Remove first image.
		 */
		Image4DType::Pointer RemoveFirstImage( const Image4DType::Pointer& input )
		{
			// reduce 4th dimension
			Image4DType::RegionType region = input->GetLargestPossibleRegion();
			Image4DType::SizeType size = region.GetSize();
			size[3]--;
			region.SetSize( size );

			// init output image
			Image4DType::Pointer output = Image4DType::New();
			output->SetRegions( region );
			output->SetOrigin( input->GetOrigin() );
			output->SetSpacing( input->GetSpacing() );
			output->Allocate();
			output->FillBuffer( 0 );

			ConstIterator4DType it( input, input->GetLargestPossibleRegion() );
			it.GoToBegin();
			it.SetDirection( 3 );

			while( ! it.IsAtEnd() )
			{
				bool firstVolume = true;

				while( ! it.IsAtEndOfLine() )
				{
					if( firstVolume )
					{
						firstVolume = false;
					}
					else
					{
						Image4DType::IndexType outputIndex = it.GetIndex();
						outputIndex[3]--;
						output->SetPixel( outputIndex, it.Get() );
					}
					++it;
				}
				it.NextLine();
			}

			return output;
		}
	};

} // end namespace dynamics

/**
 * Dynamic T1 compartment models (Logan, Patlak, Vargan, Logan_GLM, etc)
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Dynamic T1 compartment model fitting." );

	dynamics::parameters args;
	args.method = 1;
	args.relaxivity = 4.39;
	args.hematocrit = 0.38;
	args.TR = 180;
	args.maskFileName = "";
	args.verbose = false;

	// *******************************************************************************

   p.AddArgument( args.sinusFileName, "sinus" )
			->AddAlias( "s" )
			->SetDescription( "Sinus Mask (3D)" )
			->SetRequired( true );

	p.AddArgument( args.inputFileName, "input" )
			->AddAlias( "i" )
			->SetDescription( "R1 input file (first volume should be pre-bolus)" )
			->SetRequired( true );

    p.AddArgument( args.maskFileName, "mask" )
			->AddAlias( "m" )
			->SetDescription( "Mask (3D)" )
			->SetRequired( false );

	p.AddArgument( args.outputVFileName, "output-V" )
			->AddAlias( "ov" )
			->SetDescription( "Output distribution space V" )
			->SetRequired( false );
	
    p.AddArgument( args.outputKFileName, "output-K" )
			->AddAlias( "ok" )
			->SetDescription( "Output transfer constant K" )
			->SetRequired( false );
    
    p.AddArgument( args.outputResidualsFileName, "output-R" )
		    ->AddAlias( "or" )
			->SetDescription( "Output residual norm" )
			->SetRequired( false );
    
    p.AddArgument( args.method, "method" )
			->SetDescription( "Method ( 1 => Patlak (1985), 2 => Ichise (2002), 3 => Average(1,2); default: 1 )" )
			->SetRequired( false );

    p.AddArgument( args.hematocrit, "hematocrit" )
    		->SetDescription( "Venous blood hematocrit (default: 0.38)" )
    		->SetRequired( false );

    p.AddArgument( args.relaxivity, "relaxivity" )
        	->SetDescription( "Gd-DTPA relaxivity at 37 C (default: 4.39 s-1 mM-1)" )
        	->SetRequired( false );

    p.AddArgument( args.TR, "repetition-time" )
        	->SetDescription( "Time between 3D images (default: 180 sec)" )
        	->SetRequired( false );

    p.AddArgument( args.verbose, "verbose" )
            ->AddAlias( "v" )
            ->SetDescription( "Verbose" )
            ->SetRequired( false );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	dynamics::GraphicalAnalysis ga;
	ga.Run( args );

	return EXIT_SUCCESS;
}


