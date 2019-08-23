#include "dtifit.h"

namespace dtifit
{

	struct DTIFit::Impl
	{
		InputImageType::Pointer Input;
		OutputImageType::Pointer Mask;

		DesignMatrixType DesignMatrix;
		AlgoEnum::Enum Algorithm;

		OutputImageType::Pointer FA;
		OutputImageType::Pointer RA;
		OutputImageType::Pointer Trace;
		OutputImageType::Pointer L1;
		OutputImageType::Pointer L2;
		OutputImageType::Pointer L3;

		OutputVectorImageType::Pointer V1;
		OutputVectorImageType::Pointer V2;
		OutputVectorImageType::Pointer V3;

		TensorImageType::Pointer Tensor;

		RGBImageType::Pointer ColorFA;
	};

	DTIFit::DTIFit() :
		m_Impl( new Impl, true )
	{
	}

	DTIFit::~DTIFit()
	{
	}

	/**
	 * Dot product vectors (http://en.wikipedia.org/wiki/Dot_product).
	 */
	DTIFit::PixelType DTIFit::Dot( const VectorType& a, const VectorType& b )
	{
		PixelType sum = 0.0;

		for ( unsigned int k = 0; k < a.size(); k++ )
			sum += a[k] * b[k];

		return sum;
	}

	/**
	 * Dot product matrices (http://en.wikipedia.org/wiki/Dot_product).
	 */
	DTIFit::MatrixType DTIFit::Dot( const MatrixType& m1, const MatrixType& m2 )
	{
		return m1 * m2;
	}

	/**
	 * Return max value.
	 */
	DTIFit::PixelType DTIFit::Max( PixelType a, PixelType b )
	{
		if ( a < b )
			return b;
		else
			return a;
	}

	/**
	 * Rho to gamma conversion.
	 */
	void DTIFit::Gamma( const VectorType& rho, VectorType& gamma )
	{
		PixelType tmp = rho[1];

		gamma[0] = rho[0];
		gamma[1] = rho[1] * rho[1];
		gamma[2] = ( rho[2] * rho[2] + rho[4] * rho[4] );
		gamma[3] = ( rho[3] * rho[3] + rho[5] * rho[5] + rho[6] * rho[6] );
		gamma[4] = ( rho[4] * tmp );
		gamma[5] = ( rho[2] * rho[5] + rho[4] * rho[6] );
		gamma[6] = ( tmp * rho[6] );
	}

	/**
	 * Dot product: W . J.
	 */
	DTIFit::MatrixType DTIFit::WJ( const VectorType& rho, const MatrixType& W, unsigned int dim )
	{
		return Dot( W, J( rho, dim ) );
	}

	/***
	 * Jacobian.
	 */
	DTIFit::MatrixType DTIFit::J( const VectorType& rho, unsigned int dim )
	{
		double r1 = rho[1];
		double r2 = rho[2];
		double r3 = rho[3];
		double r4 = rho[4];
		double r5 = rho[5];
		double r6 = rho[6];

		MatrixType J( dim, dim, 0 );

		J[0][0] = 1.0;

		J[1][1] = ( 2.0 * r1 );
		J[4][1] = r4;
		J[6][1] = r6;

		J[2][2] = ( 2.0 * r2 );
		J[5][2] = r5;

		J[3][3] = ( 2.0 * r3 );

		J[2][4] = ( 2.0 * r4 );
		J[4][4] = r1;
		J[5][4] = r6;

		J[3][5] = ( 2.0 * r5 );
		J[5][5] = r2;

		J[3][6] = ( 2.0 * r6 );
		J[5][6] = r4;
		J[6][6] = r1;

		return J;
	}

	/**
	 * Hessian matrix.
	 */
	DTIFit::MatrixType DTIFit::Hessian( const VectorType& sihat, const VectorType& r, unsigned int dim, unsigned int N, const MatrixType& W )
	{
		MatrixType hessian( dim, dim, 0 );

		for ( unsigned int i = 0; i < N; i++ )
		{
			hessian[1][1] += -2.0 * W[i][1] * sihat[i] * r[i];
			hessian[2][2] += -2.0 * W[i][2] * sihat[i] * r[i];
			hessian[3][3] += -2.0 * W[i][3] * sihat[i] * r[i];
			hessian[1][4] += -1.0 * W[i][4] * sihat[i] * r[i];
			hessian[1][6] += -1.0 * W[i][6] * sihat[i] * r[i];
			hessian[2][5] += -1.0 * W[i][5] * sihat[i] * r[i];
		}

		hessian[4][1] = hessian[1][4];
		hessian[6][1] = hessian[1][6];
		hessian[5][2] = hessian[2][5];
		hessian[4][4] = hessian[2][2];
		hessian[5][5] = hessian[3][3];
		hessian[6][6] = hessian[3][3];
		hessian[4][6] = hessian[2][5];
		hessian[6][4] = hessian[2][5];

		return hessian;
	}

	/**
	 * Get value.
	 */
	DTIFit::PixelType DTIFit::ValueAt( const MatrixType& W, const VectorType& si, unsigned int dim, unsigned int N, const VectorType& rho )
	{
		double r = 0.0;
		double f = 0.0;
		double val1 = 0.0;

		VectorType gamma( rho.size(), 0 );
		Gamma( rho, gamma );

		for ( unsigned int i = 0; i < N; i++ )
		{
			val1 = 0.0;
			for ( unsigned int j = 0; j < dim; j++ )
				val1 += W[i][j] * gamma[j];

			r = si[i] - std::exp( val1 );
			f += 0.5 * r * r;
		}

		return f;
	}

	/**
	 * Update step.
	 */
	DTIFit::PixelType DTIFit::Update( unsigned int dim, unsigned int N, const MatrixType& W, const VectorType& si, const VectorType& rho,
			PixelType lambda, VectorType& negGrad, MatrixType& hessian_p )
	{
		MatrixType wj = WJ( rho, W, dim );

		VectorType sihat( N, 0 );
		VectorType r( N, 0 );
		double f = 0.0;
		double val1 = 0.0;
		double val2 = 0.0;
		VectorType diag( N, 0 );
		VectorType gamma( dim, 0 );

		Gamma( rho, gamma );

		for ( unsigned int i = 0; i < N; i++ )
		{
			val1 = 0.0;
			for ( unsigned int j = 0; j < dim; j++ )
				val1 += W[i][j] * gamma[j];

			sihat[i] = std::exp( val1 );
			r[i] = si[i] - sihat[i];
			f += 0.5 * r[i] * r[i];
			val2 = sihat[i] * r[i];
			diag[i] = sihat[i] * sihat[i] - val2;

			for ( unsigned int k = 0; k < dim; k++ )
				negGrad[k] += wj[i][k] * val2;
		}

		for ( unsigned int i = 0; i < dim; i++ )
		{
			for ( unsigned int j = i; j < dim; j++ )
			{
				val1 = 0.0;

				for ( unsigned int k = 0; k < N; k++ )
					val1 += wj[k][i] * diag[k] * wj[k][j];

				hessian_p[i][j] = val1;
				hessian_p[j][i] = val1;
			}

			hessian_p[i][i] += lambda;
		}

		hessian_p = hessian_p + Hessian( sihat, r, dim, N, W );

		return f;
	}

	/**
	 * Constrained weighted linear least squares solution.
	 */
	DTIFit::VectorType DTIFit::CWLLS( const MatrixType& W, const VectorType& si, bool gammaOutput )
	{
		VectorType gamma = WLLS( W, si );

		return CWLLS( gamma, gammaOutput );
	}

	/**
	 * Constrained weighted linear least squares solution.
	 */
	DTIFit::VectorType DTIFit::CWLLS( const VectorType& gamma, bool gammaOutput )
	{
		MatrixType svdTensor( 3, 3, 0 );

		svdTensor( 0, 0 ) = gamma( 1 );
		svdTensor( 0, 1 ) = gamma( 4 );
		svdTensor( 0, 2 ) = gamma( 6 );

		svdTensor( 1, 0 ) = gamma( 4 );
		svdTensor( 1, 1 ) = gamma( 2 );
		svdTensor( 1, 2 ) = gamma( 5 );

		svdTensor( 2, 0 ) = gamma( 6 );
		svdTensor( 2, 1 ) = gamma( 5 );
		svdTensor( 2, 2 ) = gamma( 3 );

		// upper triangle ...
		MatrixType R = ModifiedCholesky( svdTensor );
		VectorType psol = GetCholeskyComposition( R );

		VectorType output( 7 );
		output( 0 ) = gamma( 0 );
		output( 1 ) = psol( 0 );
		output( 2 ) = psol( 1 );
		output( 3 ) = psol( 2 );
		output( 4 ) = psol( 3 );
		output( 5 ) = psol( 4 );
		output( 6 ) = psol( 5 );

		if ( gammaOutput )
		{
			VectorType tmp( 7, 0 );
			Gamma( output, tmp );
			return tmp;
		} else
		{
			return output;
		}
	}

	/**
	 * Get upper triangular elements of input matrix as a vector.
	 */
	DTIFit::VectorType DTIFit::GetCholeskyComposition( const MatrixType& m )
	{
		unsigned int dimX = m.rows();
		unsigned int dimY = m.columns();
		unsigned int indexY = 0;
		unsigned int index = 0;
		unsigned int len = ( dimX * ( dimX + 1 ) ) / 2;

		VectorType vec( len, 0 );

		for ( unsigned int i = 0; i < dimX; dimX-- )
		{
			unsigned int px = i;
			for ( unsigned int j = indexY; j < dimY; j++ )
			{
				vec[index++] = m[px][j];
				px++;
			}

			indexY++;
		}

		return vec;
	}

	/**
	 * Constrained weighted nonlinear least squares.
	 */
	DTIFit::VectorType DTIFit::CWNLLS( const DTIFit::DesignMatrixType& W, const VectorType si )
	{
		VectorType typx( 7 );
		typx( 0 ) = 1.0;
		typx( 1 ) = 0.001;
		typx( 2 ) = 0.001;
		typx( 3 ) = 0.001;
		typx( 4 ) = 0.001;
		typx( 5 ) = 0.001;
		typx( 6 ) = 0.001;

		unsigned int dim = W.columns(); // 7
		unsigned int N = W.rows(); // 69

		VectorType x = CWLLS( W, si, false );
		VectorType x_p = x;

		double lambda = 0.0;
		VectorType negGrad( 7, 0 );
		MatrixType hessian( dim, dim, 0 );
		double f = 20.0;
		double f_p = 10.0;
		VectorType delta( dim, 0 );

		for ( int count = 0; count < 100; count++ ) // 100 iterations ...
		{
			if ( std::abs< double >( f_p - f ) < 3.5000000000000002E-013 )
			{
				double maxvalue2 = 0.0;
				double temp2 = 0.0;

				for ( int i = 0; i < 7; i++ )
				{
					temp2 = std::abs< double >( x_p[i] - x[i] ) / Max( std::abs< double >( x_p[i] ), typx[i] );
					if ( temp2 > maxvalue2 )
						maxvalue2 = temp2;
				}

				if ( maxvalue2 < 2.3099725541661634E-014 )
				{
					VectorType gamma( 7, 0 );
					Gamma( x, gamma );
					return gamma;
				}
			}

			x = x_p;

			f = Update( dim, N, W, si, x, lambda, negGrad, hessian );

			vnl_qr< double > qr( hessian );

			// hack to prevent rank deficiency warnings...
			std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original cout sbuf
			std::streambuf* cerr_sbuf = std::cerr.rdbuf(); // save original cerr sbuf
			std::ofstream fout( "/dev/null" );
			std::ofstream ferr( "/dev/null" );
			std::cout.rdbuf( fout.rdbuf() ); // redirect 'cout' to a 'fout'
			std::cerr.rdbuf( ferr.rdbuf() ); // redirect 'cerr' to a 'ferr'

			delta = qr.solve( negGrad );

			std::cout.rdbuf( cout_sbuf ); // restore the original stream buffer
			std::cerr.rdbuf( cerr_sbuf ); // restore the original stream buffer

			double alpha = 1.0;
			double fx = f;

			VectorType gradfx( dim, 0 );
			VectorType tempx( dim, 0 );

			for ( unsigned int i = 0; i < dim; i++ )
				gradfx[i] = -1 * negGrad[i];

			for ( unsigned int i = 0; i < dim; i++ )
				tempx[i] = x[i] + delta[i] * alpha;

			double fxp = ValueAt( W, si, dim, N, tempx );
			double EXTRA = 0.089999999999999997 * Dot( gradfx, delta );
			double term = fx + alpha * EXTRA;

			unsigned int INDEX = 0;

			for ( INDEX = 0; fxp > term && INDEX < 20; INDEX++ )
			{
				alpha = 0.01 * alpha;
				for ( unsigned int i = 0; i < dim; i++ )
					tempx[i] = x[i] + delta[i] * alpha;

				fxp = ValueAt( W, si, dim, N, tempx );
				term = fx + alpha * EXTRA;
			}

			if ( INDEX == 20 || fxp > term )
			{
				if ( lambda == 0.0 )
					lambda = 9.9999999999999995E-007;
				else
					lambda = 5;
				continue;
			}
			lambda = 0.0;
			for ( unsigned int i = 0; i < dim; i++ )
				x_p[i] = tempx[i];

			f_p = fxp;
		}

		// convert back to gamma ...
		VectorType gamma( 7, 0 );
		Gamma( x, gamma );

		return gamma;
	}

	/**
	 * Return absolute matrix.
	 */
	DTIFit::MatrixType DTIFit::Abs( const MatrixType& mat )
	{
		MatrixType temp = mat;

		for ( unsigned int i = 0; i < temp.rows(); i++ )
			for ( unsigned int j = 0; j < temp.columns(); j++ )
				temp[i][j] = std::abs< double >( temp[i][j] );

		return temp;
	}

	/**
	 * Return absolute vector.
	 */
	DTIFit::VectorType DTIFit::Abs( const VectorType& V )
	{
		VectorType temp = V;

		for ( unsigned int i = 0; i < temp.size(); i++ )
			temp[i] = std::abs< double >( temp[i] );

		return temp;
	}

	/**
	 * Remove diagonal values.
	 */
	DTIFit::MatrixType DTIFit::RemoveDiagonals( const MatrixType& mat )
	{
		if ( mat.rows() != mat.columns() )
		{
			std::cerr << "*** ERROR ***: could not remove diagonals, not a square matrix!" << std::endl;
			exit( EXIT_FAILURE );
		}

		MatrixType temp = mat;

		for ( unsigned int i = 0; i < temp.rows(); i++ )
			temp[i][i] = 0.0;

		return temp;
	}

	/**
	 * Return diagonals.
	 */
	DTIFit::VectorType DTIFit::GetDiagonals( const MatrixType& mat )
	{
		if ( mat.rows() != mat.columns() )
		{
			std::cerr << "*** ERROR ***: could not get diagonals, not a square matrix!" << std::endl;
			exit( EXIT_FAILURE );
		}

		VectorType temp( mat.rows(), 0 );

		for ( unsigned int i = 0; i < mat.rows(); i++ )
			temp[i] = mat[i][i];

		return temp;
	}

	/**
	 * Modified cholesky decomposition.
	 *
	 * References:
	 *
	 * (1) Gill P, Murray W, Wright MH. Practical optimization. Academic Press. 1981.
	 * (2) Nocedal J, Wright SJ. Numerical optimization. New York: Springer; 1999.
	 */
	DTIFit::MatrixType DTIFit::ModifiedCholesky( const MatrixType& A )
	{
		MatrixType a = A;

		if ( a.rows() != a.columns() )
		{
			std::cerr << "*** ERROR ***: input modified cholesky not a square matrix" << std::endl;
			exit( EXIT_FAILURE );
		}

		if ( !IsSymmetric( a ) )
		{
			std::cerr << "*** ERROR ***: input modified cholesky not a symmetric matrix" << std::endl;
			exit( EXIT_FAILURE );
		}

		double eps = MachineEpsilon();
		VectorType array( 3, 0 );
		unsigned int n = a.rows();
		VectorType diag = GetDiagonals( a );
		double gamma = ( Abs( diag ) ).max_value();
		double xi = ( Abs( RemoveDiagonals( a ) ) ).max_value();
		double delta = eps * Max( gamma + xi, 1.0 );
		array[0] = gamma;
		array[1] = xi / sqrt( pow( n, 2 ) - 1.0 );
		array[2] = eps;
		double beta = sqrt( array.max_value() );

		VectorType d( n, 0 );
		MatrixType L = GetIdentity( n );

		VectorType ccol( n, 0 );
		double cjj = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double theta = 0.0;

		for ( unsigned int j = 0; j < n; j++ )
		{
			sum = 0.0;

			for ( unsigned int k = 0; k < j; k++ )
			{
				sum += L[j][k] * d[k] * L[j][k];
			}

			cjj = a[j][j] - sum;

			if ( j < n - 1 )
			{
				for ( unsigned int i = j + 1; i < n; i++ )
				{
					temp = 0.0;

					for ( unsigned int k = 0; k < j; k++ )
					{
						temp += L[i][k] * d[k] * L[j][k];
					}

					ccol[i] = a[i][j] - temp;
				}

				theta = ( Abs( ccol ) ).max_value();

				array[0] = std::abs< double >( cjj );
				array[1] = pow( theta / beta, 2 );
				array[2] = delta;
				d[j] = ( array ).max_value();

				for ( unsigned int i = j + 1; i < n; i++ )
				{
					L[i][j] = ccol[i] / d[j];
				}
			} else
			{
				d[j] = Max( std::abs< double >( cjj ), delta );
			}
		}

		MatrixType tempL( n, n, 0 );

		for ( unsigned int i = 0; i < n; i++ )
		{
			for ( unsigned int j = 0; j < n; j++ )
			{
				tempL[i][j] = L[i][j] * sqrt( d[j] );
			}
		}

		return tempL.transpose();
	}

	/**
	 * Get identity matrix
	 */
	DTIFit::MatrixType DTIFit::GetIdentity( unsigned int n )
	{
		MatrixType I( n, n, 0 );

		for ( unsigned int r = 0; r < I.rows(); r++ )
			I( r, r ) = 1;

		return I;
	}

	/**
	 * Machine epsilon gives an upper bound
	 * on the relative error due to rounding
	 * in floating point arithmetic.
	 *
	 * http://en.wikipedia.org/wiki/Machine_epsilon
	 */
	double DTIFit::MachineEpsilon()
	{
		double tmp1 = 1.0;
		double tmp2;
		double eps;

		do
		{
			eps = tmp1;
			tmp1 /= 2;
			tmp2 = 1.0 + tmp1;
		} while ( tmp2 > 1.0 );

		return eps;
	}

	/**
	 * Check symmetry matrix.
	 */
	bool DTIFit::IsSymmetric( const MatrixType& mat )
	{
		unsigned int m = mat.rows();
		for ( unsigned int i = 0; i < m; i++ )
		{
			for ( unsigned int j = i; j < m; j++ )
				if ( mat[i][j] != mat[j][i] )
					return false;

		}
		return true;
	}

	/**
	 * Linear least square fit (checked).
	 *
	 * si is expected to be log(DWI signal)!
	 *
	 * Returns vector with:ln alpha, Dxx, Dyy, Dzz, Dxy, Dyz, Dxz.
	 */
	DTIFit::VectorType DTIFit::WLLS( const DesignMatrixType& W, const VectorType si )
	{
		if ( W.rows() != si.size() )
		{
			std::cout << "*** ERROR ***: signal array does match design matrix size!" << std::endl;
			exit( EXIT_FAILURE );
		}

		int ro = W.rows();
		int co = W.cols();

		MatrixType M( ro, co );

		// sigma matrix ...
		for ( int i = 0; i < ro; i++ )
			for ( int j = 0; j < co; j++ )
				M[i][j] = si[i] * W[i][j];

		// decompose ...
		vnl_qr< double > qr( M );

		// fit against log( DWI signal ) ...
		VectorType siLog( si.size() );

		for ( unsigned int i = 0; i < si.size(); i++ )
			siLog( i ) = log( si( i ) );

		// solve ...
		VectorType gamma = qr.solve( Multiply( si, siLog ) );

		return gamma;
	}

	/**
	 * Build full B-matrix.
	 *
	 * http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
	 */
	DTIFit::DesignMatrixType DTIFit::GetDesignMatrixVersion1( const MatrixType& r, const VectorType& b )
	{
		DesignMatrixType A( r.rows(), 7, 0 ); // 52 x 7

		for ( unsigned int i = 0; i < r.rows(); i++ )
		{
			A( i, 0 ) = 1; // ln alpha ...

			A( i, 1 ) = -b( i ) * r( i, 0 ) * r( i, 0 ); // Gx * Gx

			A( i, 2 ) = 2 * -b( i ) * r( i, 0 ) * r( i, 1 ); // 2 * Gx * Gy

			A( i, 3 ) = 2 * -b( i ) * r( i, 0 ) * r( i, 2 ); // 2 * Gx * Gz

			A( i, 4 ) = -b( i ) * r( i, 1 ) * r( i, 1 ); // Gy * Gy

			A( i, 5 ) = 2 * -b( i ) * r( i, 1 ) * r( i, 2 ); // 2 * Gy * Gz

			A( i, 6 ) = -b( i ) * r( i, 2 ) * r( i, 2 ); // Gz * Gz
		}

		return A;
	}

	/**
	 * Build full B-matrix.
	 *
	 *C.G. Koay et al / Journal of Magnetic Resonance 128 (2006) pag. 116
	 */
	DTIFit::DesignMatrixType DTIFit::GetDesignMatrixVersion2( const MatrixType& g, const VectorType& B )
	{
		DesignMatrixType A( g.rows(), 7, 0 ); // 52 x 7

		for ( unsigned int i = 0; i < g.rows(); i++ )
		{
			double b = B( i );
			double Gx = g( i, 0 );
			double Gy = g( i, 1 );
			double Gz = g( i, 2 );

			A( i, 0 ) = 1; // ln alpha ...

			A( i, 1 ) = -b * Gx * Gx;

			A( i, 2 ) = -b * Gy * Gy;

			A( i, 3 ) = -b * Gz * Gz;

			A( i, 4 ) = -2 * b * Gx * Gy;

			A( i, 5 ) = -2 * b * Gy * Gz;

			A( i, 6 ) = -2 * b * Gx * Gz;
		}

		return A;
	}

	/**
	 * Set 4D input image.
	 */
	void DTIFit::SetImage( InputImageType::Pointer input )
	{
		m_Impl->Input = input;
	}

	void DTIFit::SetMask( OutputImageType::Pointer mask )
	{
		m_Impl->Mask = mask;
	}

	/**
	 * Set fitting algorithm.
	 */
	void DTIFit::SetAlgorithm( AlgoEnum::Enum algorithm )
	{
		m_Impl->Algorithm = algorithm;
	}

	/**
	 * Set design matrix.
	 */
	void DTIFit::SetDesignMatrix( const DesignMatrixType& designMatrix )
	{
		m_Impl->DesignMatrix = designMatrix;
	}

	/**
	 * Multiply vectors and return diagonal.
	 */
	DTIFit::VectorType DTIFit::Multiply( const VectorType& a, const VectorType& b )
	{
		unsigned int N = a.size();

		MatrixType A( a.size(), 1 );
		MatrixType B( b.size(), 1 );

		for ( unsigned int i = 0; i < N; i++ )
			A( i, 0 ) = a( i );

		for ( unsigned int i = 0; i < N; i++ )
			B( i, 0 ) = b( i );

		A = B * A.transpose();

		VectorType output( N );
		for ( unsigned int i = 0; i < N; i++ )
			output( i ) = A( i, i );

		return output;
	}

	/**
	 * Return empty tensor image with input info (spacing, origin, etc).
	 */
	DTIFit::TensorImageType::Pointer DTIFit::EmptyTensorImage( OutputImageType::ConstPointer input )
	{
		TensorImageType::SpacingType spacing = input->GetSpacing();
		TensorImageType::PointType origin = input->GetOrigin();

		TensorImageType::SizeType size = input->GetLargestPossibleRegion().GetSize();
		TensorImageType::IndexType index = input->GetLargestPossibleRegion().GetIndex();
		TensorImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( index );

		TensorImageType::Pointer output = TensorImageType::New();
		output->SetSpacing( spacing );
		output->SetOrigin( origin );
		output->SetRegions( region );
		output->Allocate();
		output->FillBuffer( 0.0 );

		return output;
	}

	/**
	 * Return 2-rand tensor.
	 *
	 * ITK tensor layout =>
	 * 					| 0  1  2  |
	 *			        | X  3  4  |
	 *       			| X  X  5  |
	 */
	DTIFit::TensorPixelType DTIFit::Tensor( const VectorType& gamma )
	{
		TensorPixelType D( 0.0 );

		D( 0, 0 ) = gamma( 1 );
		D( 0, 1 ) = gamma( 4 );
		D( 0, 2 ) = gamma( 6 );
		D( 1, 1 ) = gamma( 2 );
		D( 1, 2 ) = gamma( 5 );
		D( 2, 2 ) = gamma( 3 );

		return D;
	}

	/**
	 * Fit for one voxel.
	 */
	DTIFit::VectorType DTIFit::FitVoxel( const VectorType& si )
	{
		VectorType gamma;

		if ( m_Impl->Algorithm == AlgoEnum::WeightedLinearLeastSquares )
		{
			gamma = WLLS( m_Impl->DesignMatrix, si );
		} else if ( m_Impl->Algorithm == AlgoEnum::ConstrainedWeightedLinearLeastSquares )
		{
			gamma = CWLLS( m_Impl->DesignMatrix, si, true );
		} else if ( m_Impl->Algorithm == AlgoEnum::ConstrainedWeightedNonLinearLeastSquares )
		{
			gamma = CWNLLS( m_Impl->DesignMatrix, si );
		} else
		{
			std::cerr << "*** ERROR ***: algorithm not known!" << std::endl;
			exit( EXIT_FAILURE );
		}
		return gamma;
	}

	/**
	 * Fit all voxels within mask.
	 */
	void DTIFit::Fit()
	{

		InputImageType::Pointer input = m_Impl->Input;
		OutputImageType::Pointer mask = m_Impl->Mask;

		// allocate tensor output ...
		TensorImageType::Pointer output = EmptyTensorImage( mask.GetPointer() );

		ConstIterator4DType it( input, input->GetLargestPossibleRegion() );
		ConstIterator3DType mit( mask, mask->GetLargestPossibleRegion() );
		IteratorTensorType oit( output, output->GetLargestPossibleRegion() );

		it.SetDirection( 3 );
		it.GoToBegin();
		mit.GoToBegin();
		oit.GoToBegin();

		unsigned int totalDWI = ( input->GetLargestPossibleRegion().GetSize() )[3];

		// for each series of DWI ...
		while ( !it.IsAtEnd(), !mit.IsAtEnd(), !oit.IsAtEnd() )
		{
			if ( mit.Get() != 0 )
			{
				VectorType si( totalDWI );

				while ( !it.IsAtEndOfLine() )
				{
					si( ( it.GetIndex() )[3] ) = it.Get();
					++it;
				}

				TensorPixelType tensor = Tensor( FitVoxel( si ) );
				oit.Set( tensor );
			}

			it.NextLine();
			++mit;
			++oit;
		}

		TensorImageType::RegionType region = output->GetLargestPossibleRegion();

		// compute eigen values

		OutputImageType::Pointer outputFA = OutputImageType::New();
		OutputImageType::Pointer outputRA = OutputImageType::New();
		OutputImageType::Pointer outputTrace = OutputImageType::New();

		OutputImageType::Pointer outputL1 = OutputImageType::New();
		OutputImageType::Pointer outputL2 = OutputImageType::New();
		OutputImageType::Pointer outputL3 = OutputImageType::New();

		RGBImageType::Pointer colorFA = RGBImageType::New();

		outputFA->CopyInformation( output );
		outputRA->CopyInformation( output );
		outputTrace->CopyInformation( output );

		outputL1->CopyInformation( output );
		outputL2->CopyInformation( output );
		outputL3->CopyInformation( output );

		colorFA->CopyInformation( output );

		outputFA->SetRegions( region );
		outputRA->SetRegions( region );
		outputTrace->SetRegions( region );

		outputL1->SetRegions( region );
		outputL2->SetRegions( region );
		outputL3->SetRegions( region );

		colorFA->SetRegions( region );

		outputFA->Allocate();
		outputRA->Allocate();
		outputTrace->Allocate();

		outputL1->Allocate();
		outputL2->Allocate();
		outputL3->Allocate();

		colorFA->Allocate();

		itk::ImageRegionConstIterator< TensorImageType > itTensor( output, region );

		itk::ImageRegionIterator< OutputImageType > itFA( outputFA, region );
		itk::ImageRegionIterator< OutputImageType > itRA( outputRA, region );
		itk::ImageRegionIterator< OutputImageType > itTrace( outputTrace, region );

		itk::ImageRegionIterator< OutputImageType > itL1( outputL1, region );
		itk::ImageRegionIterator< OutputImageType > itL2( outputL2, region );
		itk::ImageRegionIterator< OutputImageType > itL3( outputL3, region );

		itk::ImageRegionIterator< RGBImageType > itColorFA( colorFA, region );

		OutputVectorImageType::Pointer outputV1 = OutputVectorImageType::New();
		OutputVectorImageType::Pointer outputV2 = OutputVectorImageType::New();
		OutputVectorImageType::Pointer outputV3 = OutputVectorImageType::New();

		TensorImageType::SizeType tensorSize = region.GetSize();
		TensorImageType::SpacingType tensorSpacing = output->GetSpacing();
		TensorImageType::PointType tensorOrigin = output->GetOrigin();
		TensorImageType::DirectionType tensorDirection = output->GetDirection();

		OutputVectorImageType::RegionType vectorRegion;
		OutputVectorImageType::SizeType vectorSize;
		OutputVectorImageType::DirectionType vectorDirection;
		OutputVectorImageType::SpacingType vectorSpacing;
		OutputVectorImageType::PointType vectorOrigin;
		OutputVectorImageType::IndexType vectorIndex;

		vectorIndex.Fill( 0 );
		vectorDirection.SetIdentity();
		for ( int i = 0; i < 3; ++i )
		{
			vectorSize[i] = tensorSize[i];
			vectorSpacing[i] = tensorSpacing[i];
			vectorOrigin[i] = tensorOrigin[i];

			for ( int j = 0; j < 3; ++j )
			{
				vectorDirection[i][j] = tensorDirection[i][j];
			}
		}

		vectorSize[3] = 3;
		vectorSpacing[3] = 1;
		vectorOrigin[3] = 0;

		vectorRegion.SetIndex( vectorIndex );
		vectorRegion.SetSize( vectorSize );

		outputV1->SetRegions( vectorRegion );
		outputV1->SetSpacing( vectorSpacing );
		outputV1->SetOrigin( vectorOrigin );
		outputV1->SetDirection( vectorDirection );
		outputV1->Allocate();

		outputV2->CopyInformation( outputV1 );
		outputV2->SetRegions( vectorRegion );
		outputV2->Allocate();

		outputV3->CopyInformation( outputV1 );
		outputV3->SetRegions( vectorRegion );
		outputV3->Allocate();

		typedef itk::ImageRegionIterator< OutputVectorImageType > VectorIterator;
		std::vector< VectorIterator > vectorIterators;

		for ( int i = 0; i < 3; ++i )
		{
			vectorSize[3] = 1;
			vectorIndex[3] = i;
			vectorRegion.SetSize( vectorSize );
			vectorRegion.SetIndex( vectorIndex );

			vectorIterators.push_back( VectorIterator( outputV3, vectorRegion ) );
			vectorIterators.push_back( VectorIterator( outputV2, vectorRegion ) );
			vectorIterators.push_back( VectorIterator( outputV1, vectorRegion ) );
		}

		for ( ; !itTensor.IsAtEnd(); ++itTensor, ++itFA, ++itRA, ++itTrace, ++itL1, ++itL2, ++itL3, ++itColorFA )
		{
			const TensorPixelType& tensor = itTensor.Get();

			EigenValuesArrayType eigenValues;
			EigenVectorsMatrixType eigenVectors;
			tensor.ComputeEigenAnalysis( eigenValues, eigenVectors );

			itFA.Set( tensor.GetFractionalAnisotropy() );
			itRA.Set( tensor.GetRelativeAnisotropy() );
			itTrace.Set( tensor.GetTrace() );

			itL1.Set( eigenValues[2] );
			itL2.Set( eigenValues[1] );
			itL3.Set( eigenValues[0] );

			RGBPixelType& rgb = itColorFA.Value();
			rgb[0] = eigenValues[2] * 255. * itFA.Value();
			rgb[1] = eigenValues[1] * 255. * itFA.Value();
			rgb[2] = eigenValues[0] * 255. * itFA.Value();

			int iterator = 0;
			for ( int i = 0; i < 3; ++i ) // for each component x, y, z
			{
				for ( int j = 0; j < 3; ++j, ++iterator ) // for each
				{
					VectorIterator& it = vectorIterators[iterator];
					it.Set( eigenVectors( j, i ) );
					++it;
				}
			}
		}

		m_Impl->FA = outputFA;
		m_Impl->RA = outputRA;
		m_Impl->Trace = outputTrace;

		m_Impl->L1 = outputL1;
		m_Impl->L2 = outputL2;
		m_Impl->L3 = outputL3;

		m_Impl->V1 = outputV1;
		m_Impl->V2 = outputV2;
		m_Impl->V3 = outputV3;

		m_Impl->Tensor = output;

		m_Impl->ColorFA = colorFA;
	}

	DTIFit::OutputImageType::Pointer DTIFit::GetFA()
	{
		return m_Impl->FA;
	}

	DTIFit::OutputImageType::Pointer DTIFit::GetRA()
	{
		return m_Impl->RA;
	}

	DTIFit::OutputImageType::Pointer DTIFit::GetTrace()
	{
		return m_Impl->Trace;
	}

	DTIFit::OutputImageType::Pointer DTIFit::GetL1()
	{
		return m_Impl->L1;
	}

	DTIFit::OutputImageType::Pointer DTIFit::GetL2()
	{
		return m_Impl->L2;
	}

	DTIFit::OutputImageType::Pointer DTIFit::GetL3()
	{
		return m_Impl->L3;
	}

	DTIFit::OutputVectorImageType::Pointer DTIFit::GetV1()
	{
		return m_Impl->V1;
	}

	DTIFit::OutputVectorImageType::Pointer DTIFit::GetV2()
	{
		return m_Impl->V2;
	}

	DTIFit::OutputVectorImageType::Pointer DTIFit::GetV3()
	{
		return m_Impl->V3;
	}

	DTIFit::TensorImageType::Pointer DTIFit::GetTensor()
	{
		return m_Impl->Tensor;
	}

	DTIFit::RGBImageType::Pointer DTIFit::GetColorFA()
	{
		return m_Impl->ColorFA;
	}

	/**
	 * Run dtifit.
	 */
	void DTIFit::Run( InputImageType::Pointer input, const std::string& maskFileName, const Procparser& pp,
			const std::string& outputFileName, bool mirrorGradientTable, unsigned int algorithm )
	{
		int numberOfRO = pp.GetSize( "dro" );
		int numberOfImages = input->GetLargestPossibleRegion().GetSize()[3];

		vnl_matrix< double > table( numberOfRO, 3 );

		for ( int i = 0; i < numberOfRO; ++i )
		{
			table( i, 0 ) = pp.GetAs< double > ( "dpe", i );
			table( i, 1 ) = pp.GetAs< double > ( "dro", i );
			table( i, 2 ) = pp.GetAs< double > ( "dsl", i );
		}

		if ( numberOfImages == 2 * numberOfRO )
		{
			vnl_matrix< double > twice( numberOfRO * 2, 3 );

			if ( !mirrorGradientTable )
			{
				for ( int i = 0; i < numberOfRO; ++i )
				{
					twice.set_row( i * 2, table.get_row( i ) );
					twice.set_row( i * 2 + 1, table.get_row( i ) * -1. );
				}
			}
            else
			{
				for ( int i = 0; i < numberOfRO; ++i )
				{
					twice.set_row( i, table.get_row( i ) );
				}

				for ( int i = 0; i < numberOfRO; ++i )
				{
					twice.set_row( i + numberOfRO, table.get_row( i ) * -1. );
				}
			}

			table = twice;
			numberOfRO = 2 * numberOfRO;
		}

        if ( numberOfImages != numberOfRO )
		{
			std::cout << "WARNING: number of acquired images (" << numberOfImages << ") != number of gradient directions (" << numberOfRO
					<< ")" << std::endl;

			if ( numberOfImages > numberOfRO )
			{
				return;
			}

			table = table.extract( numberOfImages, 3 );
			numberOfRO = numberOfImages;
		}

		std::string orient = pp.GetAs< std::string > ( "orient" );

		if ( orient == "cor" )
		{
			vnl_matrix< double > permute( 3, 3 );
			permute.set_identity();

			permute( 0, 0 ) = -1;
			permute( 0, 1 ) = 0;
			permute( 0, 2 ) = 0;

			permute( 1, 0 ) = 0;
			permute( 1, 1 ) = 0;
			permute( 1, 2 ) = 1;

			permute( 2, 0 ) = 0;
			permute( 2, 1 ) = 1;
			permute( 2, 2 ) = 0;

			table = ( permute * table.transpose() ).transpose();
		} else if ( orient == "trans90" )
		{
			vnl_matrix< double > permute( 3, 3 );
			permute.set_identity();

			permute( 0, 0 ) = 0;
			permute( 0, 1 ) = 1;
			permute( 0, 2 ) = 0;

			permute( 1, 0 ) = -1;
			permute( 1, 1 ) = 0;
			permute( 1, 2 ) = 0;

			permute( 2, 0 ) = 0;
			permute( 2, 1 ) = 0;
			permute( 2, 2 ) = 1;

			table = ( permute * table.transpose() ).transpose();
		} else if ( orient == "trans" )
		{
			vnl_matrix< double > permute( 3, 3 );
			permute.set_identity();

			permute( 0, 0 ) = 1;
			permute( 0, 1 ) = 0;
			permute( 0, 2 ) = 0;

			permute( 1, 0 ) = 0;
			permute( 1, 1 ) = -1;
			permute( 1, 2 ) = 0;

			permute( 2, 0 ) = 0;
			permute( 2, 1 ) = 0;
			permute( 2, 2 ) = 1;

			table = ( permute * table.transpose() ).transpose();
		} else
		{
			std::cout << "WARNING: unsupported orientation: " << orient << std::endl;
		}

		double gdiff = pp.GetAs< double > ( "gdiff" );
		double tdelta = pp.GetAs< double > ( "tdelta" );
		double tDELTA = pp.GetAs< double > ( "tDELTA" );

		// Price and Kuchel, J Magnetic Resonance 94 (1991) 133-139 eq.5 page 137!
		double bValue = ( 2. / vnl_math::pi ) * ( 2. / vnl_math::pi ) * ( tdelta * tdelta ) * ( gdiff * gdiff ) * ( tDELTA - tdelta / 4. ) // 4 not 3 !!
				* ( 26751.98775 * 26751.98775 * 0.01 );

		std::ofstream bvals( ( outputFileName + std::string( "_bvals" ) ).c_str() );
		std::ofstream bvecs( ( outputFileName + std::string( "_bvecs" ) ).c_str() );

		VectorType b( numberOfRO );
		MatrixType r( numberOfRO, 3 );

		for ( int i = 0; i < numberOfRO; ++i )
		{
			vnl_vector< double > vector = table.get_row( i );
			b( i ) = ( vector( 0 ) * vector( 0 ) + vector( 1 ) * vector( 1 ) + vector( 2 ) * vector( 2 ) ) * bValue;
			r( i, 0 ) = vector( 0 );
			r( i, 1 ) = vector( 1 );
			r( i, 2 ) = vector( 2 );

			bvals << b( i ) << std::endl;
			bvecs << r( i, 0 ) << " " << r( i, 1 ) << " " << r( i, 2 ) << std::endl;
		}

		this->SetDesignMatrix( GetDesignMatrixVersion2( r, b ) );
		this->SetImage( input );
		this->SetMask( GetM( maskFileName, input.GetPointer() ) );

		if ( algorithm == 1 )
			this->SetAlgorithm( AlgoEnum::WeightedLinearLeastSquares );
		else if ( algorithm == 2 )
			this->SetAlgorithm( AlgoEnum::ConstrainedWeightedLinearLeastSquares );
		else if ( algorithm == 3 )
			this->SetAlgorithm( AlgoEnum::ConstrainedWeightedNonLinearLeastSquares );

		this->Fit();

		// Writer
		WriterType::Pointer writer = WriterType::New();
		VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
		TensorWriterType::Pointer tensorWriter = TensorWriterType::New();

		// FA
		writer->SetFileName( ( outputFileName + std::string( "_FA.nii.gz" ) ).c_str() );
		writer->SetInput( this->GetFA() );
		writer->Update();

		// RA
		writer->SetFileName( ( outputFileName + std::string( "_RA.nii.gz" ) ).c_str() );
		writer->SetInput( this->GetRA() );
		writer->Update();

		// Trace
		writer->SetFileName( ( outputFileName + std::string( "_Trace.nii.gz" ) ).c_str() );
		writer->SetInput( this->GetTrace() );
		writer->Update();

		// L1
		writer->SetFileName( ( outputFileName + std::string( "_L1.nii.gz" ) ).c_str() );
		writer->SetInput( this->GetL1() );
		writer->Update();

		// L2
		writer->SetFileName( ( outputFileName + std::string( "_L2.nii.gz" ) ).c_str() );
		writer->SetInput( this->GetL2() );
		writer->Update();

		// L3
		writer->SetFileName( ( outputFileName + std::string( "_L3.nii.gz" ) ).c_str() );
		writer->SetInput( this->GetL3() );
		writer->Update();

		// V1
		vectorWriter->SetFileName( ( outputFileName + std::string( "_V1.nii.gz" ) ).c_str() );
		vectorWriter->SetInput( this->GetV1() );
		vectorWriter->Update();

		// V2
		vectorWriter->SetFileName( ( outputFileName + std::string( "_V2.nii.gz" ) ).c_str() );
		vectorWriter->SetInput( this->GetV2() );
		vectorWriter->Update();

		// V3
		vectorWriter->SetFileName( ( outputFileName + std::string( "_V3.nii.gz" ) ).c_str() );
		vectorWriter->SetInput( this->GetV3() );
		vectorWriter->Update();
	}

	/**
	 * Return mask if file given, else create empty mask from input file.
	 */
	DTIFit::OutputImageType::Pointer DTIFit::GetM( const std::string& maskFileName, DTIFit::InputImageType::ConstPointer input )
	{
		if ( !maskFileName.empty() )
		{
			ImageReaderType::Pointer reader = ImageReaderType::New();
			reader->SetFileName( maskFileName );
			reader->Update();
			return reader->GetOutput();
		} else
		{
			OutputImageType::Pointer output = OutputImageType::New();
			OutputImageType::RegionType region3D;
			OutputImageType::IndexType index3D;
			OutputImageType::SizeType size3D;
			OutputImageType::SpacingType spacing3D;
			OutputImageType::PointType origin3D;

			InputImageType::RegionType region = input->GetLargestPossibleRegion();
			InputImageType::SizeType size = region.GetSize();
			InputImageType::IndexType index = region.GetIndex();
			InputImageType::SpacingType spacing = input->GetSpacing();
			InputImageType::PointType origin = input->GetOrigin();

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

			return output;
		}
	}
} // end namespace dtifit


