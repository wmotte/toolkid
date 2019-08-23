//=== File Prolog ================================================================
//	This code was developed by NASA, Goddard Space Flight Center, Code 588
//      and computer science department of University of Maryland at College Park.
//--- Contents -------------------------------------------------------------------
//	File 'IsoClus.cc' contains the implementaion of ISOCLUS algorithm in C++
//      This version can handle two implementations of the algorithm: one is the
//      straightforward algorithm (Std) and the other which differs from the original
//      algorithm only in distance norm being used and the condition used for
//      splitting (Hyb). Hyb version uses squared Euclidean distances and in step
//      8 of the algorithm uses square root of average of squared distances.
//
//     These versions are referred to as 'Std' and 'Hyb' in our paper titled "A Fast
//     Implementation of the ISODATA Clustering Algorithm" submitted to International
//     Journal of Computational Geometry and Applications.
//
//--- Description ----------------------------------------------------------------
// ISOCLUS Algorithm: This implementation is based on ISOCLUS as described in
// http://www.pcigeomatics.com/cgi-bin/pcihlp/ISOCLUS|ALGORITHM
//
// IsoClus is a modified version of ISODATA algorithm.
// For more information about the ISODATA algorithm please refer to:
//
// Tou, Julius T. and Rafael C. Gonzalez.  1974.  Pattern
// Recognition Principles.  Addison-Wesley Publishing Co.
//
//-- Notes:-----------------------------------------------------------------------
// This program was developed as part of a cooperation between Code 588, and
// Code 935 (Dr. Jacqueline LeMoigne), and computer science department of
// University of Maryland (Dr. David Mount and Dr. Nathan Netanyahu). Future plans
// are to improve this algorithm in the similar to the way that Dr. Mount improved
// k-means algorithm by using special data structures (ie: kd-tree, ANN library elements)
// and methods. For  information about the efficient implementation of the k-means
// algorithm done by UMCP please refer to:
// http://www.cs.umd.edu/~mount/pubs.html
// T. Kanungo, D. M. Mount, N. Netanyahu, C. Piatko, R. Silverman, and A. Y. Wu
// , ``An Efficient k-Means Clustering Algorithm: Analysis and Implementation''
//  IEEE Trans.PAMI, 24 (2002), 881-892.
//
//-- Development History:--------------------------------------------------------
//   Date	      Author		    Reference
//   Description
//
//   May 24, 2002     Nargess Memarsadeghi  NASA GSFC, Code 588
//   Initial implementation
//
//--- DISCLAIMER---------------------------------------------------------------
//	This software is provided "as is" without any warranty of any kind, either
//	express, implied, or statutory, including, but not limited to, any
//	warranty that the software will conform to specification, any implied
//	warranties of merchantability, fitness for a particular purpose, and
//	freedom from infringement, and any warranty that the documentation will
//	conform to the program, or any warranty that the software will be error
//	free.
//
//	In no event shall NASA be liable for any damages, including, but not
//	limited to direct, indirect, special or consequential damages, arising out
//	of, resulting from, or in any way connected with this software, whether or
//	not based upon warranty, contract, tort or otherwise, whether or not
//	injury was sustained by persons or property or otherwise, and whether or
//	not loss was sustained from or arose out of the results of, or use of,
//	their software or services provided hereunder.
//
//--- Warning-------------------------------------------------------------------
//    This software is property of the National Aeronautics and Space
//    Administration.  Unauthorized use or duplication of this software is
//    strictly prohibited.  Authorized users are subject to the following
//    restrictions:
//    *   Neither the author, their corporation, nor NASA is responsible for
//        any consequence of the use of this software.
//    *   The origin of this software must not be misrepresented either by
//        explicit claim or by omission.
//    *   Altered versions of this software must be plainly marked as such.
//    *   This notice may not be removed or altered.
//
//=== End File Prolog=============================================================

#include <iostream.h>
#include <fstream.h>
#include <string>
#include <ctype.h>
#include <time.h>

#include <vector>
#include "Image.h"
#include "Point.h"
#include "KM_ANN.h"
#include "KMrand.h"

// Overview:
// ---------
// The program is run as follows:
//
//      IsoClus < input > output
//
// where the 'input' file contains a list of directives as described
// below.  Directives consist of a directive name, followed by list of
// arguments (depending on the directive).  Arguments and directives are
// separated by white space (blank, tab, and newline).  String arguments
// are not quoted, and consist of a string of nonwhite chacters.  A
// character "#" denotes a comment.  The following characters up to
// the end of line are ignored.  Comments may only be inserted between
// directives (not within the argument list of a directive).
//
// The program in most cases assumes that use put valid data for
// input parameters and avoids detailed error checkings.
//
// Directives:
// ----------
// NUMCLUS <int>
// NSAM <int>
// NUMBANDS <int>
// SAMPRM <int>
// STDV <double>
// LUMP <double>
// SQUARED <int>
// MAXPAIR <int>
// MAXITER <int>
// rows <int>
// cols <int>
// input_images <string*>
// output <string>
//------------------------------------------------------------------------
//  Default execution parameters
//------------------------------------------------------------------------
const int DEF_num_clus = 16; // number of desired clusters
const int DEF_NSAM = 262144; // number of samples to collect for iterative clustering.
const int DEF_num_bands = 0; // Number of the image bands (dimension of the space). This is the
// required input from user and this defaut value is used and termination
// condition if the user does not provide any value for it.
const int DEF_num_samples = 5; // minimum nubmer of samples in a cluster
const int DEF_max_pair = 5; // maximum number of pairs of cluster centers that can be lumped.
const int DEF_max_iter = 20; // maximum number of iterations to run the clustering algorithm.


const int DEF_rows = 0; // Number of image rows.  Required input.  Default=0 used
// for error checking if not provided.
const int DEF_cols = 0; // Number of image columns.  Required input.
// Default=0 used for error checking if not provided.
const float DEF_stdv = 10; // standard deviation within a cluster.
const float DEF_data_std = 0.005;
const float DEF_lump = 1; // the lumping parameter
const int DEF_squared = 0; // An integer to determine whether or not to calculate
// squared distances in step 5 (as opposed to regular
// distances).

string* DEF_input_images = NULL; // name of the input image files.
string DEF_output = "classified";// the name of the file in which the classified image will be stored.
// checking if not provided.

const int DEF_seed = 0; //seed for random generation of points
// (when generating synthetic data)
const int DEF_sample_seed = 0; // seed for randomly sampling centers and points
ostream* const DEF_out = &cout; //Standard output stream.
ostream* const DEF_err = &cerr; //error output stream.
istream* const DEF_in = &cin; //Standard input stream.

//-----------------------------------------------------------------------
// Gloabal constants.
//-----------------------------------------------------------------------
enum Isoerr {
	Isowarn = 0, Isoabort = 1
}; // what to do in case of error
//------------------------------------------------------------------------
//  IsoClus global variables - Execution options
//------------------------------------------------------------------------
int NUMCLUS;
int NSAM;
int NUMBANDS;
int SAMPRM;
int MAXPAIR;
int MAXITER;
int rows;
int cols;
double STDV;
double data_std;
double LUMP;
int SQUARED;
string* input_images;
string output;
string distribution = "clus_gauss"; // if user selects to generate data, this can be either
// of the following: uniform, gauss, laplace, co_gauss,
// co_laplace, clus_gauss, clus_orth_flats,clus_ellipsoids,
// multi_clus
int seed; //for random number generator (when generating synthetic data).
int sample_seed; // for random sampling of centers and points
bool is_rand;

istream* IsoIn; //standard input stream
ostream* IsoOut; //standard output stream
ostream* IsoErr;

ifstream inStream; // input stream
ofstream outStream; // output stream

//------------------------------------------------------------------------
//  Initialize global parameters
//------------------------------------------------------------------------
static void initGlobals() {

	NUMCLUS = DEF_num_clus;
	NSAM = DEF_NSAM;
	NUMBANDS = DEF_num_bands;
	SAMPRM = DEF_num_samples;
	MAXPAIR = DEF_max_pair;
	MAXITER = DEF_max_iter;
	rows = DEF_rows;
	cols = DEF_cols;
	STDV = DEF_stdv;
	data_std = DEF_data_std;
	LUMP = DEF_lump;
	SQUARED = DEF_squared;
	input_images = DEF_input_images;
	output = DEF_output;

	IsoIn = DEF_in;
	IsoOut = DEF_out;
	IsoErr = DEF_err;
	kmIdum = -DEF_seed;
	sample_seed = DEF_sample_seed;
	is_rand = false;

}
//------------------------------------------------------------------------
//  IsoExit - exit from program
//      This is used mostly because Borland C++ deletes the command
//      prompt window immediately on exit.  This keeps until
//      the user verifies deletion.
//------------------------------------------------------------------------
void IsoExit( int status ) // exit program
{
#ifdef _BORLAND_C
	char ch;
	if (IsoIn == &cin) { // input from std in
		cerr << "Hit return to continue..." << endl;
		IsoIn->get(ch);
	}
#endif
	exit( status );
}
//------------------------------------------------------------------------
//  IsoError - print error message
//      If Isoerr is Isoabort we also abort the program.
//------------------------------------------------------------------------

void IsoError( const string &msg, Isoerr level ) {
	if ( level == Isoabort ) {
		*IsoErr << "IsoClus: ERROR------->" << msg << "<-------------ERROR" << endl;
		*IsoOut << "IsoClus: ERROR------->" << msg << "<-------------ERROR" << endl;
		IsoExit( 1 );
	} else {
		*IsoErr << "IsoClus: WARNING----->" << msg << "<-------------WARNING" << endl;
		*IsoOut << "IsoClus: WARNING----->" << msg << "<-------------WARNING" << endl;
	}
}
//------------------------------------------------------------------------
// getCmdArgs - get and process command line arguments
//
//      Syntax:
//      kmltest [-i infile] [-o outfile]
//
//      where:
//          infile              directive input file
//          outfile             output file
//
//      If file is not specified, then the standard input and standard
//      output are the defaults.
//------------------------------------------------------------------------
void getCmdArgs( int argc, char *argv[] ) {
	int i = 1;
	while ( i < argc ) { // read arguments
		if ( !strcmp( argv[i], "-i" ) ) { // -i option
			inStream.open( argv[++i], ios::in ); // open input file
			if ( !inStream ) {
				IsoError( "Cannot open input file", Isoabort );
			}
			IsoIn = &inStream; // make this input stream
		} else if ( !strcmp( argv[i], "-o" ) ) { // -o option
			outStream.open( argv[++i], ios::out );// open output file
			if ( !outStream ) {
				IsoError( "Cannot open output file", Isoabort );

			}
			IsoOut = &outStream; // make this output stream
		} else { // illegal syntax
			*IsoErr << "Syntax:\n" << "kmltest [-i infile] [-o outfile]\n" << "    where:\n" << "    infile         directive input file\n"
					<< "    outfile        output file\n";
			exit( 1 ); // exit
		}
		i++;
	}
}

//--------------------------------------------------------------------------
// Note: This function is written by Dr. Mount amd was originally part of his
//       kmltest.cpp program.
//--------------------------------------------------------------------------
static bool skipComment( istream &in ) {
	char ch = 0;
	// skip whitespace
	do {
		in.get( ch );
	} while ( isspace( ch ) && !in.eof() );
	while ( ch == '#' && !in.eof() ) { // comment?
		// skip to end of line
		do {
			in.get( ch );
		} while ( ch != '\n' && !in.eof() );
		// skip whitespace
		do {
			in.get( ch );
		} while ( isspace( ch ) && !in.eof() );
	}
	if ( in.eof() )
		return false; // end of file
	in.putback( ch ); // put character back
	return true;
}
//------------------------------------------------------------------------
// Note: This function is written by Dr. Mount amd was originally part of his
//       kmltest.cpp program.
// getDirective - skip comments and read next directive
//      Returns true if directive read, and false if eof seen.
//------------------------------------------------------------------------

static bool getDirective( istream &in, string &directive ) {
	if ( !skipComment( in ) ) // skip comments
		return false; // found eof along the way?
	in >> directive; // read directive
	return true;
}
//----------------------------------------------------------------------
// elapsedTime
// Utility for computing elapsed time.
//----------------------------------------------------------------------

#ifndef CLOCKS_PER_SEC                  // define clocks-per-second if needed
#define CLOCKS_PER_SEC          1000000 // (should be in time.h)
#endif

inline double elapsedTime( clock_t start ) {
	return double( clock() - start ) / double( CLOCKS_PER_SEC);
}
//--------------------------------------------------------------------------
int main( int argc, char** argv ) {

	string directive; //input directive
	string strArg; // string argument

	KMpointArray all = NULL;

	initGlobals();
	getCmdArgs( argc, argv );
	IsoOut->precision( 4 ); // output precision

	// read input directive
	while ( getDirective( *IsoIn, directive ) ) {
		//----------------------------------------------------------------
		//  Read options
		//----------------------------------------------------------------
		if ( directive == "NUMCLUS" ) {
			*IsoIn >> NUMCLUS;

		} else if ( directive == "NSAM" ) {
			*IsoIn >> NSAM;
		} else if ( directive == "NUMBANDS" ) {
			*IsoIn >> NUMBANDS;
		} else if ( directive == "SAMPRM" ) {
			*IsoIn >> SAMPRM;
		} else if ( directive == "STDV" ) {
			*IsoIn >> STDV;
		} else if ( directive == "data_std" ) {

			*IsoIn >> data_std;
		} else if ( directive == "LUMP" ) {
			*IsoIn >> LUMP;
		} else if ( directive == "SQUARED" ) {
			*IsoIn >> SQUARED;
		} else if ( directive == "MAXPAIR" ) {
			*IsoIn >> MAXPAIR;
		} else if ( directive == "MAXITER" ) {
			*IsoIn >> MAXITER;
		} else if ( directive == "rows" ) {
			*IsoIn >> rows;
		} else if ( directive == "cols" ) {
			*IsoIn >> cols;
		} else if ( directive == "input_images" ) {
			//the user must have specified number of bands in the image before the names of the input data.
			if ( NUMBANDS == 0 ) {
				IsoError( "You should specify NUMBANDS before input_images.", Isoabort );
			}
			input_images = new string[NUMBANDS];
			skipComment( *IsoIn );
			for ( int i = 0; i < NUMBANDS; i++ ) {
				*IsoIn >> input_images[i];
			}
		} else if ( directive == "output" ) {
			*IsoIn >> output;
		} else if ( directive == "seed" ) {
			*IsoIn >> kmIdum;
			kmIdum = -kmIdum;

		} else if ( directive == "sample_seed" ) {
			*IsoIn >> sample_seed;
		} else if ( directive == "gen_data_pts" ) {

			all = kmAllocPts( rows * cols, NUMBANDS );

			// double sep=1.0;
			// kmClusGaussPts(all, rows*cols, NUMBANDS, NUMCLUS,true,data_std, sep);
			double std_m = 1.0 / pow( NUMCLUS, 1 / (double) NUMBANDS );
			double std_std = std_m / 3;
			STDV = 2 * std_m;
			kmClusGaussPts2( all, rows * cols, NUMBANDS, NUMCLUS, true, std_m, std_std );

			is_rand = true;

		} else {
			*IsoErr << "Directive: " << directive << "\n";
			*IsoErr << "Unknown directive";
			exit( 1 );
		}

	}
	// if the user have not specified any of the required inputs, generate
	// an error message and exit the program.

	if ( NUMBANDS == 0 ) {
		IsoError( "You need to specify the directive called NUMBANDS, that is number of bands in the image", Isoabort );

	} else if ( rows == 0 ) {
		IsoError( "You need to specify the directive called rows, that is number of rows in the image.", Isoabort );
	} else if ( cols == 0 ) {
		IsoError( "You need to specify the directive called cols, that is number of columns in the image.", Isoabort );

	}

	//if the user have not specified any value for NSAM, and the default
	//value is greater than the image size, set it to image size;

	if ( NSAM > rows * cols )
		NSAM = rows * cols;

	//------------------------------------------------------------
	int iter = 0;
	double exec_time = 0;

	Image IMG = Image( rows, cols, NUMBANDS, NUMCLUS, SAMPRM );

	if ( is_rand ) {

		IMG.setPoints( all );

		// STEP 1:
		//read image
	} else {
		*IsoOut << "\n\tReading image files." << endl;

		IMG.readImage( input_images );
	}

	*IsoOut << "\tSelecting initial centers by sampling." << endl;
	IMG.sampleCenters();

	*IsoOut << "\tSampling " << NSAM << " image points for iterative clustering." << endl;
	IMG.samplePoints( NSAM );

	clock_t start = clock(); // start the clock

	for ( iter = 1; iter <= MAXITER; iter++ ) {
		*IsoOut << "\nIteration Number " << iter << " :" << endl;
		if ( iter == MAXITER ) {
			*IsoOut << "\tPerform the last iterative clustering on all points.\n";
			IMG.preFinalClustering();
		}
		do {
			//STEP 2:

			//calculate the distances of all points from each center to look up
			// these distances later

			IMG.CalculateDistances();

			//now, for each point we need to figure out to which cluster it belongs to.
			//basically we just need to find the minimum number in each column of the distances
			*IsoOut << "\tPut points into clusters." << endl;
			IMG.PutInCluster();

			//STEP 3:
			IMG.PostAnalyzeClusters();

			//STEP 4:
			*IsoOut << "\tUpdate centers by calculating the average point in each cluster." << endl;
			IMG.UpdateCenters();

		} while ( IMG.WasDeleted() );

		//need to update distances since in the last iteration centers have modified.

		IMG.CalculateDistances();

		IMG.PutInCluster();
		//STEP 5:
		IMG.CalculateAverageDistances();

		//STEP 6:

		IMG.OverallAverageDistances();

		//STEP 7:
		int next_step = 8;
		if ( iter == MAXITER ) {
			LUMP = 0;
			next_step = 11;
		}

		else if ( IMG.getNumCenters() <= ( NUMCLUS / 2 ) ) {
			next_step = 8;
		} else if ( ( iter % 2 == 0 ) || ( IMG.getNumCenters() >= 2 * NUMCLUS ) ) {
			next_step = 11;
		}

		switch ( next_step ) {
		case 8: {
			//STEP 8:

			IMG.CalculateSTDVector();

			//STEP 9:

			IMG.CalculateVmax();

			//STEP 10:
			// the vector to_split will contain integers that represent the cluster numbers
			// that need to be split.

			vector< int > to_split = IMG.ShouldSplit( STDV );

			if ( to_split.size() != 0 ) {
				IMG.Split( to_split );
				//we need to substract one if it was the last iteration because otherwise we
				//we will exit the loop without updating clusters.
				if ( iter == MAXITER )
					iter = iter - 1;
				break; //go to step 2
			}

		} //CASE 8

		case 11: {
			//STEP 11:

			IMG.ComputeCenterDistances();

			//STEP 12:
			vector< PairDistanceNode > to_lump = IMG.FindLumpCandidates( LUMP, MAXPAIR );

			//STEP 13:
			if ( to_lump.size() != 0 )
				IMG.Lump( to_lump );

		} //CASE 11


		}// SWITCH


	}// FOR LOOP


	exec_time = elapsedTime( start ); // get elapsed time

	if ( !is_rand )
		IMG.writeClassifiedImage( output );

	IMG.generateReport();
	*IsoOut << "Algorithm's run time: " << exec_time << " CPU seconds." << endl;

	*IsoOut << "total overall dist " << IMG.getDistortions() << endl;
	if ( is_rand )
		kmDeallocPts( all );

}

