#include "dkifitCommon.h"

namespace dki
{
	/**
	 * Constructor.
	 */
	DkiFit::DkiFit( const parameters& args )
	{
		InitData( args.inputFileName );
		InitMask( args.maskFileName );

		m_sep = args.sep;
		ReadBVals( args.bvalsFileName );
		ReadBVecs( args.bvecsFileName );
		RemainDWIOnly();

		InitA_D();
		InitA_K();
		InitDesignMatrix();
	}

	/**
	 * Init data (read input, allocate output).
	 */
	void DkiFit::InitData( const std::string& inputFileName )
	{
		if ( !inputFileName.empty() )
		{
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName( inputFileName );
			reader->Update();
			m_Input = reader->GetOutput();

			OutputImageType::RegionType region = m_Input->GetLargestPossibleRegion();
			OutputImageType::SizeType size = region.GetSize();
			OutputImageType::IndexType index = region.GetIndex();
			OutputImageType::SpacingType spacing = m_Input->GetSpacing();
			OutputImageType::PointType origin = m_Input->GetOrigin();

			size[3] = m_OUTPUT_IMAGES;
			index[3] = 0;

			region.SetSize( size );
			region.SetIndex( index );

			// create
			m_Output = OutputImageType::New();

			// set
			m_Output->SetRegions( region );
			m_Output->SetSpacing( spacing );
			m_Output->SetOrigin( origin );
			m_Output->Allocate();
			m_Output->FillBuffer( 0 );

		} else
		{
			std::cerr << "Could not read input!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	/**
	 * Init mask if file given, else create empty mask from input file.
	 */
	void DkiFit::InitMask( const std::string& maskFileName )
	{
		if ( !maskFileName.empty() )
		{
			MaskReaderType::Pointer reader = MaskReaderType::New();
			reader->SetFileName( maskFileName );
			reader->Update();
			m_Mask = reader->GetOutput();
		} else
		{
			MaskImageType::Pointer output = MaskImageType::New();
			MaskImageType::RegionType region3D;
			MaskImageType::IndexType index3D;
			MaskImageType::SizeType size3D;
			MaskImageType::SpacingType spacing3D;
			MaskImageType::PointType origin3D;

			ImageType::RegionType region = m_Input->GetLargestPossibleRegion();
			ImageType::SizeType size = region.GetSize();
			ImageType::IndexType index = region.GetIndex();
			ImageType::SpacingType spacing = m_Input->GetSpacing();
			ImageType::PointType origin = m_Input->GetOrigin();

			size3D[0] = size[0];
			size3D[1] = size[1];
			size3D[2] = size[2];
			index3D[0] = index[0];
			index3D[1] = index[1];
			index3D[2] = index[2];
			origin3D[0] = origin[0];
			origin3D[1] = origin[1];
			origin3D[2] = origin[2];
			spacing3D[0] = spacing[0];
			spacing3D[1] = spacing[1];
			spacing3D[2] = spacing[2];

			region3D.SetSize( size3D );
			region3D.SetIndex( index3D );

			// set
			output->SetRegions( region3D );
			output->SetSpacing( spacing3D );
			output->SetOrigin( origin3D );
			output->Allocate();
			output->FillBuffer( 1 );

			m_Mask = output;
		}
	}

	/**
	 * Write output images.
	 */
	void DkiFit::Write( const std::string& outputFileName )
	{
		// Writer
		WriterType::Pointer writer = WriterType::New();

		// DTI tensor (7 elements)
		writer->SetFileName( outputFileName.c_str() );
		writer->SetInput( m_Output );
		writer->Update();
	}

	/**
	 * Read bvals from file.
	 */
	void DkiFit::ReadBVals( const std::string& bvalsFileName )
	{
		// init
		VectorType b( GetNumberOfRows( bvalsFileName ), 0 );

		// open data ...
		std::ifstream in( bvalsFileName.c_str() );
		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read from: " << bvalsFileName << std::endl;
			exit( EXIT_FAILURE );
		}

		std::string line;
		unsigned int rowIndex = 0;

		while ( getline( in, line ) )
		{
			TokType tok( line, boost::char_separator< char >( m_sep.c_str() ) );
			try
			{
				for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
				{
					b( rowIndex ) = boost::lexical_cast< PixelType >( *id );
				}
				rowIndex++;
			} catch ( boost::bad_lexical_cast& e )
			{
				std::cout << "*** WARNING ***: could not parse " << bvalsFileName << std::endl;
				e.what();
				exit( EXIT_FAILURE );
			}
		}

		in.close();

		m_bvals = b;
	}

	void DkiFit::ReadBVecs( const std::string& bvecsFileName )
	{
		// init
		MatrixType r( GetNumberOfRows( bvecsFileName ), 3, 0 );

		// open data ...
		std::ifstream in( bvecsFileName.c_str() );
		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read from: " << bvecsFileName << std::endl;
			exit( EXIT_FAILURE );
		}

		std::string line;
		unsigned int rowIndex = 0;

		while ( getline( in, line ) )
		{
			TokType tok( line, boost::char_separator< char >( m_sep.c_str() ) );
			try
			{
				unsigned int colIndex = 0;
				for ( TokType::iterator id = tok.begin(); id != tok.end(); ++id )
				{
					r( rowIndex, colIndex ) = boost::lexical_cast< PixelType >( *id );
					colIndex++;
				}
				rowIndex++;
			} catch ( boost::bad_lexical_cast& e )
			{
				std::cout << "*** WARNING ***: could not parse " << bvecsFileName << std::endl;
				std::cout << e.what() << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		in.close();

		m_bvecs = r;
	}

	/**
	 * Return number of training data points.
	 */
	unsigned int DkiFit::GetNumberOfRows( const std::string& inputFile )
	{
		std::ifstream in( inputFile.c_str() );

		if ( in.fail() )
		{
			std::cerr << "*** ERROR ***: could not read from: " << inputFile << std::endl;
			exit( EXIT_FAILURE );
		}
		std::string line;
		unsigned int result = 0;

		while ( getline( in, line ) )
			result++;

		in.close();
		return result;
	}

	/**
	 * Init A_D matrix (is equal to B-matrix without b-value multiplication).
	 *
	 * Elements:
	 * ---------
	 * [Dxx, Dyy, Dzz, Dxy, Dyz, Dxz]
	 */
	void DkiFit::InitA_D()
	{
		// check ...
		if ( m_bvals.size() != m_bvecs.rows() )
		{
			std::cerr << "*** ERROR ***: bvals and bvecs size does not match!" << std::endl;
			exit( EXIT_FAILURE );
		}

		m_A_D = MatrixType( m_M, 6, 0 ); // M x 6

		for ( unsigned int i = 0; i < m_M; i++ )
		{
			PixelType Gx = m_bvecs( i, 0 );
			PixelType Gy = m_bvecs( i, 1 );
			PixelType Gz = m_bvecs( i, 2 );

			m_A_D( i, 0 ) = Gx * Gx;
			m_A_D( i, 1 ) = Gy * Gy;
			m_A_D( i, 2 ) = 2 * Gx * Gy;
			m_A_D( i, 3 ) = Gz * Gz;
			m_A_D( i, 4 ) = 2 * Gy * Gz;
			m_A_D( i, 5 ) = 2 * Gx * Gz;
		}
	}

	/**
	 * Init A_K matrix.
	 *
	 * Elements:
	 * ---------
	 * [W_1111, W_2222, W_3333,
	 * 	W_1112, W_2223, W_2221,
	 *  W_2223, W_3331, W_3332,
	 *  W_1122, W_1133, W_2233,
	 * 	W_1123, W_1223, W_1233]
	 */
	void DkiFit::InitA_K()
	{
		// check ...
		if ( m_bvals.size() != m_bvecs.rows() )
		{
			std::cerr << "*** ERROR ***: bvals and bvecs size does not match!" << std::endl;
			exit( EXIT_FAILURE );
		}

		m_A_K = MatrixType( m_M, 15, 0 ); // M x 15

		for ( unsigned int i = 0; i < m_M; i++ )
		{
			PixelType Gx = m_bvecs( i, 0 );
			PixelType Gy = m_bvecs( i, 1 );
			PixelType Gz = m_bvecs( i, 2 );

			m_A_K( i, 0 ) = Gx * Gx * Gx * Gx; // W_1111
			m_A_K( i, 1 ) = Gy * Gy * Gy * Gy; // W_2222
			m_A_K( i, 2 ) = Gz * Gz * Gz * Gz; // W_3333

			m_A_K( i, 3 ) = 4 * Gx * Gx * Gx * Gy; // W_1112
			m_A_K( i, 4 ) = 4 * Gx * Gx * Gx * Gz; // W_2223
			m_A_K( i, 5 ) = 4 * Gx * Gy * Gy * Gy; // W_2221
			m_A_K( i, 6 ) = 4 * Gy * Gy * Gy * Gz; // W_2223
			m_A_K( i, 7 ) = 4 * Gx * Gz * Gz * Gz; // W_3331
			m_A_K( i, 8 ) = 4 * Gy * Gz * Gz * Gz; // W_3332

			m_A_K( i, 9 ) = 6 * Gx * Gx * Gy * Gy; // W_1122
			m_A_K( i, 10 ) = 6 * Gx * Gx * Gz * Gz; // W_1133
			m_A_K( i, 11 ) = 6 * Gy * Gy * Gz * Gz; // W_2233

			m_A_K( i, 12 ) = 12 * Gx * Gx * Gy * Gz; // W_1123
			m_A_K( i, 13 ) = 12 * Gx * Gy * Gy * Gz; // W_1223
			m_A_K( i, 14 ) = 12 * Gx * Gy * Gz * Gz; // W_1233
		}
	}

	/**
	 * Init Design matrix A in: AX = b.
	 */
	void DkiFit::InitDesignMatrix()
	{
		m_DesignMatrix = MatrixType( m_A_D.rows(), m_A_D.cols() + m_A_K.cols() ); // M x 21

		// -b * A_D

		for( unsigned int r = 0; r < m_M; r++)
			for( unsigned int c = 0; c < m_A_D.cols(); c++ )
				m_DesignMatrix( r, c ) = - m_bvals( r ) * m_A_D( r, c );

		// 1/6 * -b^2 A_K

		for( unsigned int r = 0; r < m_M; r++)
			for( unsigned int c = 0; c < m_A_K.cols(); c++ )
				m_DesignMatrix( r, c + 6 ) = ( 1. / 6. ) * m_bvals( r ) * m_bvals( r ) * m_A_K( r, c );
	}

	/**
	 * 1. Remove entries without b-values == 0 from bvals and bvecs.
	 * 2. Set M to total DWIs.
	 * 3. Return vector with DWIs indices.
	 */
	void DkiFit::RemainDWIOnly()
	{
		IndexVectorType indices( m_bvals.size(), 1 );
		m_M = m_bvals.size();

		// first get indices (and m_M)

		for ( unsigned int i = 0; i < m_bvals.size(); i++ )
		{
			if ( m_bvals( i ) == 0 )
				indices( i ) = 0;
		}

		m_M = indices.sum();

		// create new matrices

		VectorType bvals = VectorType( m_M );
		MatrixType bvecs = MatrixType( m_M, 3 );

		unsigned int ind = 0;

		for ( unsigned i = 0; i < m_bvals.size(); i++ )
		{
			if ( m_bvals( i ) != 0 )
			{
				bvals( ind ) = m_bvals( i );
				bvecs.set_row( ind, m_bvecs.get_row( i ) );
				++ind;
			}
		}

		m_bvals = bvals;
		m_bvecs = bvecs;
		m_nonzero_indices = indices;
	}

	/**
	 * 1. Average S_0.
	 * 2. Remove S_0s.
	 * 3. Return ln(S_N / S_0)
	 */
	VectorType DkiFit::LogDWIOnlyData( const VectorType data )
	{
		PixelType S0 = 0;
		unsigned int N = 0;

		 // average S_0

		for( unsigned int i = 0; i < data.size(); i++ )
		{
			if( m_nonzero_indices( i ) == 0 )
			{
				S0 += data( i );
				N++;
			}
		}

		S0 /= N;

		// ln( S_N / S_0 )

		VectorType output( m_M, 0 );

		unsigned int ind = 0;

		for( unsigned int i = 0; i < data.size(); i++ )
			if( m_nonzero_indices( i ) != 0 )
			{
				output( ind ) = vcl_log( data( i ) / S0 );
				ind++;
			}

		return output;
	}

	/**
	 * Return median value.
	 */
	PixelType DkiFit::Median( const VectorType& V )
	{
	    VectorType vals = V;
		VectorType::iterator medianIter = vals.begin() + vals.size() / 2;
	    vcl_nth_element( vals.begin(), medianIter, vals.end() );

	    return *medianIter;
	}

	/**
	 * If one row contains only 0 elements, the determinant will be 0 as well.
	 * In that case reset the diagnonal to 0.1
	 */
	MatrixType DkiFit::CheckDeterminant( const MatrixType& XWX )
	{
		if ( std::abs( vnl_determinant< PixelType >( XWX ) ) <= 1e-8 )
		{
			return XWX;
		}
		else
		{
			MatrixType ret = XWX;
			for ( unsigned int i = 0; i < ret.rows(); i++ )
				ret( i, i ) = ret( i, i ) + 0.1;
			return ret;
		}
	}

	/**
	 * Resturn residuals for weighted fit.
	 */
	VectorType DkiFit::Residuals( const MatrixType& A, const DiagMatrixType& W, const VectorType& b, const VectorType& x )
	{
		return W * ( b - A * x ); // TODO check with Umesh!
	}

	/**
	 * Resturn residuals.
	 */
	VectorType DkiFit::Residuals( const MatrixType& A, const VectorType& b, const VectorType& x )
	{
		return b - A * x; // Ax = b
	}

	/**
	 * Private Constructor.
	 */
	DkiFit::DkiFit()
	{
	}

} // end namespace dki
