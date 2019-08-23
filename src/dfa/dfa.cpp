//#include <iomanip>
//#include <cmath>
//#include <ctime>
//#include <algorithm>
//#include <vector>
//#include <iterator>
//#include <functional>
//#include <numeric>
#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SWAP(a,b) {temp = (a); (a) = (b); (b) = temp;}

// *******************************************************
extern "C"
{
	/* Function prototypes. */
	void start( double *data, int numberOfPoints, float lower, float higher );
	int input( void );
	int rscale( int minbox, int maxbox, double boxratio );
	void dfa( double *seq, int npts, int nfit, int *rs, int nr, int sw );
	void setup( void );
	void cleanup( void );
	void help( void );
	double polyfit( double **x, double *y, int ndat, int nfit );
	void error( char error_text[] );
	double *vector( int nl, int nh );
	int *ivector( int nl, int nh );
	int *lvector( int nl, int nh );
	double **matrix( int nrl, int nrh, int ncl, int nch );
	void free_vector( double *v, int nl, int nh );
	void free_ivector( int *v, int nl, int nh );
	void free_lvector( int *v, int nl, int nh );
	void free_matrix( double **m, int nrl, int nrh, int ncl, int nch );

	/* Global variables. */
	char *pname; /* this program's name (for use in error messages) */
	double *seq; /* input data buffer; allocated and filled by input() */
	int *rs; /* box size array; allocated and filled by rscale() */
	double *mse; /* fluctuation array; allocated by setup(), filled by dfa() */
	int iflag = 1; /* integrate the input data if non-zero */
	int nfit = 2; /* order of the regression fit, plus 1 */
	int nr; /* number of box sizes */

	/**
	 * START.
	 */
	void start( double *data, int numberOfPoints, float lower, float higher )
	{

		int i;
		int sw = 0;
		int minbox = 0;
		int maxbox = 0;
		int npts;
		int temp;

		/* Read and interpret the command line. */
		pname = const_cast< char* > ( "dfa" );

		nfit = 2; // d ...

		minbox = lower;
		maxbox = higher;
		sw = 0; // sliding-window

		/* Allocate and fill the input data array seq[]. */
		//npts = input();
		seq = data;
		npts = numberOfPoints;

		/* Set minimum and maximum box sizes. */
		if ( minbox < 2 * nfit )
		minbox = 2 * nfit;
		if ( maxbox == 0 || maxbox> npts / 4 )
		maxbox = npts / 4;
		if ( minbox> maxbox )
		{
			SWAP(minbox, maxbox);
			if ( minbox < 2 * nfit )
			minbox = 2 * nfit;
		}

		/* Allocate and fill the box size array rs[].  rscale's third argument
		 specifies that the ratio between successive box sizes is 2^(1/8). */
		nr = rscale( minbox, maxbox, pow( 2.0, 1.0 / 8.0 ) );

		/* Allocate memory for dfa() and the functions it calls. */
		setup();

		/* Measure the fluctuations of the detrended input data at each box size
		 using the DFA algorithm; fill mse[] with these results. */
		dfa( seq, npts, nfit, rs, nr, sw );

		/* Output the results. */
		for ( i = 1; i <= nr; i++ )
		printf( "%g %g\n", log10( (double) rs[i] ), log10( mse[i] ) / 2.0 );

		/* Release allocated memory. */
		cleanup();
		exit( 0 );
	}

	/* Read input data, allocating and filling seq[], integrating if iflag != 0.
	 Following the convention used for other arrays in this program, seq[0] is
	 unused, and the first point is stored in seq[1].  The return value is the
	 number of points read.

	 This function allows the input buffer to grow as large as necessary, up to
	 the available memory (assuming that a int int is large enough to address
	 any memory location).  Note that the integration is done using double
	 precision arithmetic to avoid complete loss of precision when the integrated
	 data reach large amplitudes.  */
	int input()
	{
		//int maxdat = 0, npts = 0;
		int maxdat = 0, npts = 0;
		double y, yp = 0.0;

		while ( scanf( "%lf", &y ) == 1 )
		{
			if ( ++npts >= maxdat )
			{
				double *s;

				maxdat += 50000; /* allow the input buffer to grow (the
				 increment is arbitrary) */
				if ( ( s = (double*) realloc( seq, maxdat * sizeof(double) ) ) == NULL )
				{
					fprintf( stderr, "%s: insufficient memory, truncating input at row %d\n", pname, npts );
					break;
				}
				seq = s;
			}
			seq[npts] = iflag ? ( yp += y ) : y;
		}

		if ( npts < 1 )
			error( "no data read" );
		return ( npts );
	}

	int rslen; /* length of rs[] */

	/* rscale() allocates and fills rs[], the array of box sizes used by dfa()
	 below.  The box sizes range from (exactly) minbox to (approximately) maxbox,
	 and are arranged in a geometric series such that the ratio between
	 consecutive box sizes is (approximately) boxratio.  The return value is
	 the number of box sizes in rs[].
	 */
	int rscale( int minbox, int maxbox, double boxratio )
	{
		int ir, n;
		int rw;

		/* Determine how many scales are needed. */
		rslen = log10( maxbox / (double) minbox ) / log10( boxratio ) + 1.5;

		/* Thanks to Peter Domitrovich for pointing out that a previous version
		 of the above calculation undercounted the number of scales in some
		 situations. */
		rs = lvector( 1, rslen );

		for ( ir = 1, n = 2, rs[1] = minbox; n <= rslen && rs[n - 1] < maxbox; ir++ )
		{
			if ( ( rw = minbox * pow( boxratio, ir ) + 0.5 ) > rs[n - 1] )
			{
				rs[n++] = rw;
			}
		}

		if ( rs[--n] > maxbox )
		{
			--n;
		}

		return ( n );
	}

	double **x; /* matrix of abscissas and their powers, for polyfit(). */

	/* Detrended fluctuation analysis
	 seq:	input data array
	 npts:	number of input points
	 nfit:	order of detrending (2: linear, 3: quadratic, etc.)
	 rs:		array of box sizes (uniformly distributed on log scale)
	 nr:		number of entries in rs[] and mse[]
	 sw:		mode (0: non-overlapping windows, 1: sliding window)
	 This function returns the mean squared fluctuations in mse[].
	 */
	void dfa( double *seq, int npts, int nfit, int *rs, int nr, int sw )
	{
		int i, boxsize, inc, j;
		double stat;

		for ( i = 1; i <= nr; i++ )
		{
			boxsize = rs[i];

			if ( sw )
			{
				inc = 1;
				stat = (int) ( npts - boxsize + 1 ) * boxsize;
			}
			else
			{
				inc = boxsize;
				stat = (int) ( npts / boxsize ) * boxsize;
			}
			for ( mse[i] = 0.0, j = 0; j <= npts - boxsize; j += inc )
			{

				double polyfitResult = polyfit( x, seq + j, boxsize, nfit );
				mse[i] += polyfitResult;
			}
			mse[i] /= stat;
		}
	}

	/* workspace for polyfit() */
	double *beta, **covar, **covar0;
	int *indxc, *indxr, *ipiv;

	/* This function allocates workspace for dfa() and polyfit(), and sets
	 x[i][j] = i**(j-1), in preparation for polyfit(). */
	void setup()
	{
		int i;
		int j, k;

		beta = vector( 1, nfit );
		covar = matrix( 1, nfit, 1, nfit );
		covar0 = matrix( 1, nfit, 1, nfit );
		indxc = ivector( 1, nfit );
		indxr = ivector( 1, nfit );
		ipiv = ivector( 1, nfit );
		mse = vector( 1, nr );
		x = matrix( 1, rs[nr], 1, nfit );
		for ( i = 1; i <= rs[nr]; i++ )
		{
			printf("i %d\n", i );
			x[i][1] = 1.0;
			x[i][2] = i;
			for ( j = 3; j <= nfit; j++ )
				x[i][j] = x[i][j - 1] * i;
		}
	}

	/* This function frees all memory previously allocated by this program. */
	void cleanup()
	{
		free_matrix( x, 1, rs[nr], 1, nfit );
		free_vector( mse, 1, nr );
		free_ivector( ipiv, 1, nfit );
		free_ivector( indxr, 1, nfit );
		free_ivector( indxc, 1, nfit );
		free_matrix( covar0, 1, nfit, 1, nfit );
		free_matrix( covar, 1, nfit, 1, nfit );
		free_vector( beta, 1, nfit );
		free_lvector( rs, 1, rslen ); /* allocated by rscale() */
		free( seq ); /* allocated by input() */
	}

	static char *help_strings[] =
	{ "usage: %s [OPTIONS ...]\n", "where OPTIONS may include:", " -d K             detrend using a polynomial of degree K",
			"                   (default: K=1 -- linear detrending)", " -h               print this usage summary",
			" -i               input series is already integrated", " -l MINBOX        smallest box width (default: 2K+2)",
			" -s               sliding window DFA", " -u MAXBOX        largest box width (default: NPTS/4)",
			"The standard input should contain one column of data in text format.",
			"The standard output is two columns: log(n) and log(F) [base 10 logarithms],",
			"where n is the box size and F is the root mean square fluctuation.", NULL};

	void help( void )
	{
		int i;

		(void) fprintf( stderr, help_strings[0], pname );
		for ( i = 1; help_strings[i] != NULL; i++ )
			(void) fprintf( stderr, "%s\n", help_strings[i] );
	}

	/* polyfit() is based on lfit() and gaussj() from Numerical Recipes in C
	 (Press, Teukolsky, Vetterling, and Flannery; Cambridge U. Press, 1992).  It
	 fits a polynomial of degree (nfit-1) to a set of boxsize points given by
	 x[1...boxsize][2] and y[1...boxsize].  The return value is the sum of the
	 squared errors (chisq) between the (x,y) pairs and the fitted polynomial.
	 */
	double polyfit( double **x, double *y, int boxsize, int nfit )
	{
		int icol, irow, j, k;
		double big, chisq, pivinv, temp;
		int i;
		static int pboxsize = 0;

		/* This block sets up the covariance matrix.  Provided that boxsize
		 never decreases (which is true in this case), covar0 can be calculated
		 incrementally from the previous value. */
		if ( pboxsize != boxsize )
		{ /* this will be false most of the time */
			if ( pboxsize > boxsize ) /* this should never happen */
				pboxsize = 0;
			if ( pboxsize == 0 ) /* this should be true the first time only */
				for ( j = 1; j <= nfit; j++ )
					for ( k = 1; k <= nfit; k++ )
						covar0[j][k] = 0.0;
			for ( i = pboxsize + 1; i <= boxsize; i++ )
				for ( j = 1; j <= nfit; j++ )
					for ( k = 1, temp = x[i][j]; k <= j; k++ )
						covar0[j][k] += temp * x[i][k];
			for ( j = 2; j <= nfit; j++ )
				for ( k = 1; k < j; k++ )
					covar0[k][j] = covar0[j][k];
			pboxsize = boxsize;
		}
		for ( j = 1; j <= nfit; j++ )
		{
			beta[j] = ipiv[j] = 0;
			for ( k = 1; k <= nfit; k++ )
				covar[j][k] = covar0[j][k];
		}
		for ( i = 1; i <= boxsize; i++ )
		{
			beta[1] += ( temp = y[i] );
			beta[2] += temp * i;
		}
		if ( nfit > 2 )
			for ( i = 1; i <= boxsize; i++ )
				for ( j = 3, temp = y[i]; j <= nfit; j++ )
					beta[j] += temp * x[i][j];
		for ( i = 1; i <= nfit; i++ )
		{
			big = 0.0;
			for ( j = 1; j <= nfit; j++ )
				if ( ipiv[j] != 1 )
					for ( k = 1; k <= nfit; k++ )
					{
						if ( ipiv[k] == 0 )
						{
							if ( ( temp = covar[j][k] ) >= big || ( temp = -temp ) >= big )
							{
								big = temp;
								irow = j;
								icol = k;
							}
						} else if ( ipiv[k] > 1 )
							error( "singular matrix" );
					}
			++( ipiv[icol] );
			if ( irow != icol )
			{
				for ( j = 1; j <= nfit; j++ )
					SWAP(covar[irow][j], covar[icol][j]);
				SWAP(beta[irow], beta[icol]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if ( covar[icol][icol] == 0.0 )
				error( "singular matrix" );
			pivinv = 1.0 / covar[icol][icol];
			covar[icol][icol] = 1.0;
			for ( j = 1; j <= nfit; j++ )
				covar[icol][j] *= pivinv;
			beta[icol] *= pivinv;
			for ( j = 1; j <= nfit; j++ )
				if ( j != icol )
				{
					temp = covar[j][icol];
					covar[j][icol] = 0.0;
					for ( k = 1; k <= nfit; k++ )
						covar[j][k] -= covar[icol][k] * temp;
					beta[j] -= beta[icol] * temp;
				}
		}
		chisq = 0.0;
		if ( nfit <= 2 )
			for ( i = 1; i <= boxsize; i++ )
			{
				temp = beta[1] + beta[2] * i - y[i];
				chisq += temp * temp;
			}
		else
			for ( i = 1; i <= boxsize; i++ )
			{
				temp = beta[1] + beta[2] * i - y[i];
				for ( j = 3; j <= nfit; j++ )
					temp += beta[j] * x[i][j];
				chisq += temp * temp;
			}
		return ( chisq );
	}

	/* The functions below are based on those of the same names in Numerical
	 Recipes (see above). */
	void error( char error_text[] )
	{
		fprintf( stderr, "%s: %s\n", pname, error_text );
		exit( 1 );
	}

	double *vector( int nl, int nh )
	/* allocate a double vector with subscript range v[nl..nh] */
	{
		double *v = (double *) malloc( ( size_t )( ( nh - nl + 2 ) * sizeof(double) ) );
		if ( v == NULL )
			error( "allocation failure in vector()" );
		return ( v - nl + 1 );
	}

	int *ivector( int nl, int nh )
	/* allocate an int vector with subscript range v[nl..nh] */
	{
		int *v = (int *) malloc( ( size_t )( ( nh - nl + 2 ) * sizeof(int) ) );
		if ( v == NULL )
			error( "allocation failure in ivector()" );
		return ( v - nl + 1 );
	}

	int *lvector( int nl, int nh )
	/* allocate a int int vector with subscript range v[nl..nh] */
	{
		int *v = (int *) malloc( ( size_t )( ( nh - nl + 2 ) * sizeof(int) ) );
		if ( v == NULL )
			error( "allocation failure in lvector()" );
		return ( v - nl + 1 );
	}

	double **matrix( int nrl, int nrh, int ncl, int nch )
	/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
	{
		int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
		double **m;

		/* allocate pointers to rows */
		m = (double **) malloc( ( size_t )( ( nrow + 1 ) * sizeof(double*) ) );
		if ( !m )
			error( "allocation failure 1 in matrix()" );
		m += 1;
		m -= nrl;

		/* allocate rows and set pointers to them */
		m[nrl] = (double *) malloc( ( size_t )( ( nrow * ncol + 1 ) * sizeof(double) ) );
		if ( !m[nrl] )
			error( "allocation failure 2 in matrix()" );
		m[nrl] += 1;
		m[nrl] -= ncl;

		for ( i = nrl + 1; i <= nrh; i++ )
			m[i] = m[i - 1] + ncol;

		/* return pointer to array of pointers to rows */
		return ( m );
	}

	void free_vector( double *v, int nl, int nh )
	/* free a double vector allocated with vector() */
	{
		free( v + nl - 1 );
	}

	void free_ivector( int *v, int nl, int nh )
	/* free an int vector allocated with ivector() */
	{
		free( v + nl - 1 );
	}

	void free_lvector( int *v, int nl, int nh )
	/* free a int int vector allocated with lvector() */
	{
		free( v + nl - 1 );
	}

	void free_matrix( double **m, int nrl, int nrh, int ncl, int nch )
	/* free a double matrix allocated by matrix() */
	{
		free( m[nrl] + ncl - 1 );
		free( m + nrl - 1 );
	}

}

// *******************************************************


/**
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC Utrecht, NL.
 * @date: 19-11-2009
 * This method was first proposed in:
 * Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL. Mosaic
 * organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.  [Available
 * on-line at http://prola.aps.org/abstract/PRE/v49/i2/p1685_1]
 *
 * A detailed description of the algorithm and its application to physiologic
 * signals can be found in:
 * Peng C-K, Havlin S, Stanley HE, Goldberger AL. Quantification of scaling
 * exponents and crossover phenomena in nonstationary heartbeat time series.
 * Chaos 1995;5:82-87. [Abstract online at http://www.ncbi.nlm.nih.gov/entrez/-
 * query.fcgi?cmd=Retrieve&db=PubMed&list_uids=11538314&dopt=Abstract]
 *
 * If you use this program in support of published research, please include a
 * citation of at least one of the two references above, as well as the standard
 * citation for PhysioNet:
 * Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
 * Mietus JE, Moody GB, Peng CK, Stanley HE.  PhysioBank, PhysioToolkit, and
 * Physionet: Components of a New Research Resource for Complex Physiologic
 * Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages;
 * http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13).
 */
class DetrendedFluctuationAnalysis
{

public:

	typedef double ValueType;
	typedef std::vector< ValueType > VectorType;

	/**
	 * Detrended fluctuation analysis.
	 */
	void run( const std::string& inputFileName, float lower, float upper )
	{
		// [ 1 ] read input from text file ...
		double* data;
		int inputSize = getInput( inputFileName, data );

		// [ 2 ] insert in C-start function ...

		// TODO set proper default values ...
		int d = 1;

		// call C-old main function ...
		start( data, inputSize, lower, upper );
	}

protected:

	/**
	 * Return input from given text file as vector.
	 */
	int getInput( const std::string& input, double* inputDataPointer )
	{
		std::ifstream inFile;

		inFile.open( input.c_str() );
		if ( !inFile )
		{
			std::cout << "*** ERROR ***: Unable to open: " << input << "." << std::endl;
			exit( EXIT_FAILURE );
		}

		double x;
		VectorType output;

		while ( inFile >> x )
		{
			output.push_back( x );
		}
		inFile.close();

		/**
		 * Negative values will be converted to complex numbers in matlab,
		 * but not with the stl ...
		 *
		 * No support is given (yet) for complex number mle.
		 */
		if ( *( std::min_element( output.begin(), output.end() ) ) < 0 )
		{
			std::cerr << "*** ERROR ***: Negative input not supported!" << std::endl;
			exit( EXIT_FAILURE );
		}

		inputDataPointer = new double[output.size()];
		std::copy( output.begin(), output.end(), inputDataPointer );

		// dereference iterator and get address....
		//VectorType::iterator it = output.begin();
		// set pointer ...
		//inputDataPointer = &*it;

		return output.size();
	}
};

// ************************************************************************************

/**
 * Throw error if required option is not specified.
 */
void required_option( const boost::program_options::variables_map& vm, const std::string& required_option )
{
	if ( vm.count( required_option ) == 0 )
		throw std::logic_error( "Option: '" + required_option + "' is required!" );
}

/**
 * Fit powerlaw to list of numbers.
 */
int main( int argc, char* argv[] )
{
	namespace po = boost::program_options;

	// application description ...
	std::string description = "Detrended fluctuation analysis.\n";

	// options ...
	std::string input;

	double lower;
	double upper;

	try
	{

		po::options_description desc( "Available options" );

		desc.add_options()

		( "input,i", po::value< std::string >( &input ), "string: input file with distribution values in column format." )

		( "lower,l", po::value< double >( &lower )->default_value( 0.0 ), "float: lower." )

		( "upper,u", po::value< double >( &upper )->default_value( 0.0 ), "float: upper." )

		( "help,h", "bool: produce help message." );

		po::variables_map vm;
		po::store( po::parse_command_line( argc, argv, desc ), vm );
		po::notify( vm );

		// help message ...
		if ( vm.count( "help" ) )
		{
			std::cout << argv[0] << ": " << description << std::endl;
			std::cout << desc << "\n";

			return EXIT_SUCCESS;
		}

		// required options ...
		required_option( vm, "input" );

		// run application ...
		DetrendedFluctuationAnalysis dfa;

		dfa.run( input, lower, upper );

	} catch ( std::exception& e )
	{
		std::cerr << "*** ERROR ***: " << e.what() << "\n";
		std::cerr << "Use \"" << argv[0] << " --help\" for information about application usage." << std::endl;

		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
