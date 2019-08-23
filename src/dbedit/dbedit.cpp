#include "tkdCmdParser.h"
#include "itkImage.h"
#include <cmath>
#include "anautils.h"
#include "zfstream.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <zlib.h>

/**
 * LUTS recognized by FSLView in header-record: hist.aux_file="<LUT>".
 */
struct LUTS
{
	static const int Greyscale = 0;
	static const int Red_Yellow = 1;
	static const int Blue_Lightblue = 2;
	static const int Red = 3;
	static const int Blue = 4;
	static const int Green = 5;
	static const int Yellow = 6;
	static const int Pink = 7;
	static const int Hot = 8;
	static const int Cool = 9;
	static const int Copper = 10;
	static const int render1 = 11;
	static const int render1t = 12;
	static const int render2 = 13;
	static const int render2t = 14;
	static const int render3 = 15;
	static const int MGH_Cortical = 16;
	static const int MGH_Subcortical = 17;
	static const int Random_Rainbow = 18;
};

/* macros for data type byte swapping */
#define ICONV(a) (*(int *)cconv((char*) &(a), sizeof(int), 1))
#define FCONV(a) (*(float *)cconv((char*) &(a), sizeof(float), 1))
#define FCONVM(a,b) ((float *)cconv((char*) (a), sizeof(float), b))
#define SICONV(a) (*(short int *)cconv((char*) &(a), sizeof(short int), 1))
#define SICONVM(a,b) ((short int *)cconv((char*) (a), sizeof(short int), b))
#define mcalloc(a,b) calloc(a,b)
#define hdrloc(_xyz_) ((int)(&(hdr._xyz_))-(int)(&hdr))
#define WORD -1
#define INT -2
#define FLOAT -3
#define CHAR -4
#define BYTE -5

/**
 * Char convert.
 */
char* cconv( char *dbuf, int elsize, size_t elnum )
{
	int dd2, dm1;
	char *k, *m, *ep, tv;
	dd2 = ( elsize + 1 ) / 2;
	dm1 = elsize - 1;
	k = dbuf;
	m = k + dm1;
	ep = dbuf + ( elnum * elsize );
	while ( k < ep )
	{
		tv = *k;
		*k++ = *m;
		*m-- = tv;
		if ( k >= m )
		{
			k += dd2;
			m = k + dm1;
		}
	}
	return dbuf;
}

/**
 * Header convert.
 */
struct dsr headconv( struct dsr hdr )
{
	hdr.hk.sizeof_hdr = ICONV(hdr.hk.sizeof_hdr);
	hdr.hk.extents = ICONV(hdr.hk.extents);
	hdr.hk.session_error = SICONV(hdr.hk.session_error);
	SICONVM(hdr.dime.dim, 8);
	hdr.dime.unused1 = SICONV(hdr.dime.unused1);
	hdr.dime.datatype = SICONV(hdr.dime.datatype);
	hdr.dime.bitpix = SICONV(hdr.dime.bitpix);
	hdr.dime.dim_un0 = SICONV(hdr.dime.dim_un0);

	FCONVM(hdr.dime.pixdim, 8);
	hdr.dime.vox_offset = FCONV(hdr.dime.vox_offset);
	hdr.dime.funused1 = FCONV(hdr.dime.funused1);
	hdr.dime.funused2 = FCONV(hdr.dime.funused2);
	hdr.dime.funused3 = FCONV(hdr.dime.funused3);
	hdr.dime.cal_max = FCONV(hdr.dime.cal_max);
	hdr.dime.cal_min = FCONV(hdr.dime.cal_min);

	hdr.dime.compressed = ICONV(hdr.dime.compressed);
	hdr.dime.verified = ICONV(hdr.dime.verified);
	hdr.dime.glmax = ICONV(hdr.dime.glmax);
	hdr.dime.glmin = ICONV(hdr.dime.glmin);

	SICONVM(hdr.hist.originator,5);

	hdr.hist.views = ICONV(hdr.hist.views);
	hdr.hist.vols_added = ICONV(hdr.hist.vols_added);
	hdr.hist.start_field = ICONV(hdr.hist.start_field);
	hdr.hist.field_skip = ICONV(hdr.hist.field_skip);
	hdr.hist.omax = ICONV(hdr.hist.omax);
	hdr.hist.omin = ICONV(hdr.hist.omin);
	hdr.hist.smax = ICONV(hdr.hist.smax);
	hdr.hist.smin = ICONV(hdr.hist.smin);

	return hdr;
}

/**
 * Header info.
 */
typedef struct headinfo
{
	char txt[64];
	int flag;
	int type;
	float val;
	char valtxt[64];
	int loc;
} Headinfo;

static struct dsr hdr;

/**
 * Header information info.
 */
static Headinfo iinfo[] =
{ "hk.sizeof_hdr", 0, INT, 0, "", 0, "hk.data_type", 0, 10, 0, "", 4, "hk.db_name", 0, 18, 0, "", 14, "hk.extents", 0, INT, 0, "", 32,
		"hk.session_error", 0, WORD, 0, "", 36, "hk.regular", 0, BYTE, 0, "", 38, "hk.hkey_un0", 0, BYTE, 0, "", 39, "dime.dim[0]", 0,
		WORD, 0, "", 40 + 0, "dime.dim[1]", 0, WORD, 0, "", 40 + 2, "dime.dim[2]", 0, WORD, 0, "", 40 + 4, "dime.dim[3]", 0, WORD, 0, "",
		40 + 6, "dime.dim[4]", 0, WORD, 0, "", 40 + 8, "dime.dim[5]", 0, WORD, 0, "", 40 + 10, "dime.dim[6]", 0, WORD, 0, "", 40 + 12,
		"dime.dim[7]", 0, WORD, 0, "", 40 + 14, "dime.vox_units", 0, 4, 0, "", 40 + 16, "dime.cal_units", 0, 8, 0, "", 40 + 20,
		"dime.unused1", 0, WORD, 0, "", 40 + 24, "dime.datatype", 0, WORD, 0, "", 40 + 30, "dime.bitpix", 0, WORD, 0, "", 40 + 32,
		"dime.dim_un0", 0, WORD, 0, "", 40 + 34, "dime.pixdim[0]", 0, FLOAT, 0, "", 40 + 36, "dime.pixdim[1]", 0, FLOAT, 0, "", 40 + 40,
		"dime.pixdim[2]", 0, FLOAT, 0, "", 40 + 44, "dime.pixdim[3]", 0, FLOAT, 0, "", 40 + 48, "dime.pixdim[4]", 0, FLOAT, 0, "", 40 + 52,
		"dime.pixdim[5]", 0, FLOAT, 0, "", 40 + 56, "dime.pixdim[6]", 0, FLOAT, 0, "", 40 + 60, "dime.pixdim[7]", 0, FLOAT, 0, "", 40 + 64,
		"dime.vox_offset", 0, FLOAT, 0, "", 40 + 68, "dime.funused1", 0, FLOAT, 0, "", 40 + 72, "dime.funused2", 0, FLOAT, 0, "", 40 + 76,
		"dime.cal_max", 0, FLOAT, 0, "", 40 + 84, "dime.cal_min", 0, FLOAT, 0, "", 40 + 88, "dime.compressed", 0, INT, 0, "", 40 + 92,
		"dime.verified", 0, INT, 0, "", 40 + 96, "dime.glmax", 0, INT, 0, "", 40 + 100, "dime.glmin", 0, INT, 0, "", 40 + 104,
		"hist.descrip", 0, 80, 0, "", 148 + 0, "hist.aux_file", 0, 24, 0, "", 148 + 80, "hist.orient", 0, BYTE, 0, "", 148 + 104,
		"hist.originator", 0, 10, 0, "", 148 + 105, "spm.origx", 0, WORD, 0, "", 148 + 105, "spm.origy", 0, WORD, 0, "", 148 + 107,
		"spm.origz", 0, WORD, 0, "", 148 + 109, "hist.generated", 0, 10, 0, "", 148 + 115, "hist.scannum", 0, 10, 0, "", 148 + 125,
		"hist.patient_id", 0, 10, 0, "", 148 + 135, "hist.exp_date", 0, 10, 0, "", 148 + 145, "hist.hist_un0", 0, 3, 0, "", 148 + 165,
		"hist.views", 0, INT, 0, "", 148 + 168, "hist.vols_added", 0, INT, 0, "", 148 + 172, "hist.start_field", 0, INT, 0, "", 148 + 176,
		"hist.field_skip", 0, INT, 0, "", 148 + 180, "hist.omax", 0, INT, 0, "", 148 + 184, "hist.omin", 0, INT, 0, "", 148 + 188,
		"hist.smax", 0, INT, 0, "", 148 + 192, "hist.smin", 0, INT, 0, "", 148 + 196 };

#define niinfo sizeof(iinfo)/sizeof(Headinfo)

/**
 * assign Headinfo fields to header struct.
 */
void assign( struct dsr* header, Headinfo iinfo )
{
	union
	{
		char b[4];
		short w;
		int i;
		float f;
	} un;
	char *hp;
	hp = (char *) header;

	if ( iinfo.type == WORD )
	{
		un.w = (short) iinfo.val;
		hp[iinfo.loc + 0] = un.b[0];
		hp[iinfo.loc + 1] = un.b[1];
	} else if ( iinfo.type == INT )
	{
		un.i = (int) iinfo.val;
		hp[iinfo.loc + 0] = un.b[0];
		hp[iinfo.loc + 1] = un.b[1];
		hp[iinfo.loc + 2] = un.b[2];
		hp[iinfo.loc + 3] = un.b[3];
	} else if ( iinfo.type == FLOAT )
	{
		un.f = (float) iinfo.val;
		hp[iinfo.loc + 0] = un.b[0];
		hp[iinfo.loc + 1] = un.b[1];
		hp[iinfo.loc + 2] = un.b[2];
		hp[iinfo.loc + 3] = un.b[3];
	} else if ( iinfo.type == CHAR || iinfo.type == BYTE )
	{
		hp[iinfo.loc] = (char) iinfo.val;
	} else if ( iinfo.type > 0 )
	{
		strncpy( &( hp[iinfo.loc] ), iinfo.valtxt, iinfo.type - 1 );
	}
}

/**
 * DBEDIT
 *
 * @author: wim@invivonmr.uu.nl; Image Sciences Institute, UMC Utrecht.
 */
class DBEdit
{

public:

	/**
	 * Run dbedit for all images.
	 */
	void Run( const std::vector< std::string >& inputFileNames, float cal_min, float cal_max, int label )
	{
		Init();

		SetMinMaxValue( cal_min, cal_max );

		SetLUT( label );

		for ( unsigned int i = 0; i < inputFileNames.size(); i++ )
		{
			ChangeHeader( inputFileNames.at( i ), cal_min, cal_max, label );
		}
	}

protected:

	/**
	 * Init header.
	 */
	void Init()
	{
		hdr.dime.dim[0] = 4;
		hdr.dime.dim[4] = 1;
		hdr.hk.regular = 'r';
		hdr.hk.sizeof_hdr = sizeof(struct dsr);
		hdr.dime.dim[0] = 4;
		hdr.dime.dim[4] = 1;
		hdr.hk.regular = 'r';
		hdr.hk.sizeof_hdr = sizeof(struct dsr);
	}

	/**
	 * Set min and max header fields.
	 */
	void SetMinMaxValue( float cal_min, float cal_max )
	{
		// if max; do nothing...
		if ( cal_max < std::numeric_limits< float >::max() )
		{
			int value = 32; // dime.cal_min
			iinfo[value].flag |= 2;
			iinfo[value].val = cal_max;
		}

		// if min; do nothing...
		if ( cal_min > std::numeric_limits< float >::min() )
		{
			int value = 33; // dime.cal_min
			iinfo[value].flag |= 2;
			iinfo[value].val = cal_min;
		}
	}

	/**
	 * Set LUT header.
	 */
	void SetLUT( int label )
	{
		// change label ...
		if ( label != 0 )
		{
			int value = 39;
			iinfo[value].flag |= 2;

			// Red-Yellow
			if( label == 0 )
				strncpy( iinfo[value].valtxt, "Greyscale", 63 );
			else if ( label == 1 )
				strncpy( iinfo[value].valtxt, "Red-Yellow", 63 );
			else if ( label == 2 )
				strncpy( iinfo[value].valtxt, "Blue-Lightblue", 63 );
			else if ( label == 3 )
				strncpy( iinfo[value].valtxt, "Red", 63 );
			else if ( label == 4 )
				strncpy( iinfo[value].valtxt, "Green", 63 );
			else if ( label == 5 )
				strncpy( iinfo[value].valtxt, "Blue", 63 );
			else if ( label == 6 )
				strncpy( iinfo[value].valtxt, "Yellow", 63 );
			else if ( label == 7 )
				strncpy( iinfo[value].valtxt, "Pink", 63 );
			else if ( label == 8 )
				strncpy( iinfo[value].valtxt, "Cool", 63 );
			else if ( label == 9 )
				strncpy( iinfo[value].valtxt, "Hot", 63 );
			else if ( label == 10 )
				strncpy( iinfo[value].valtxt, "render1", 63 );
			else if ( label == 11 )
				strncpy( iinfo[value].valtxt, "render2", 63 );
			else if ( label == 12 )
				strncpy( iinfo[value].valtxt, "render3", 63 );
			else if ( label == 13 )
				strncpy( iinfo[value].valtxt, "MGH Cortical", 63 );
			else if ( label == 14 )
				strncpy( iinfo[value].valtxt, "MGH Subcortical", 63 );
			else if ( label == 15 )
				strncpy( iinfo[value].valtxt, "Random-Rainbow", 63 );
			else
			{
				std::cout << "Label not recognized! (Field cleaned)." << std::endl;
				strncpy( iinfo[value].valtxt, "", 63 );
			}
		}
	}

	/**
	 * Return true if gzipped.
	 */
	bool IsGzip( const std::string& inputFileName )
	{
		if ( inputFileName.rfind( ".nii.gz" ) != std::string::npos )
		{
			return true;
		} else
		{
			return false;
		}
	}

	/**
	 * Change cal_min, cal_max or label field in header.
	 */
	void ChangeHeader( const std::string& inputFileName, float cal_min, float cal_max, int label )
	{
		bool hswf = false; // header swap flag.
		bool zipped = IsGzip( inputFileName );

		FILE *hdrfp;
		gzFile fdi, fdo;

		if ( zipped )
		{
			fdi = gzopen( inputFileName.c_str(), "rb" );
			if ( !fdi )
			{
				std::cerr << "failed to open " << inputFileName << std::endl;
				exit( EXIT_FAILURE );
			}
		} else
		{
			hdrfp = fopen( inputFileName.c_str(), "r+" );
			if ( !hdrfp )
			{
				std::cerr << "failed to open " << inputFileName << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		std::vector< unsigned char > data;
		if ( zipped )
		{
			gzread( fdi, (char *) ( &hdr ), sizeof(struct dsr) );

			unsigned char buffer[1024];
			int n;
			while ( ( n = gzread( fdi, buffer, sizeof(unsigned char) * 1024 ) ) != 0 )
			{
				for ( int i = 0; i < n; ++i )
				{
					data.push_back( buffer[i] );
				}
			}
			gzclose( fdi );
		} else
		{
			// open header ...
			if ( fread( (char *) ( &hdr ), sizeof(struct dsr), 1, hdrfp ) != 1 )
			{
				std::cerr << "Warning: header " << inputFileName << " to small" << std::endl;
			}
		}

		/* check for swapped header */
		if ( ( hdr.hk.sizeof_hdr == 1543569408 ) || ( hdr.dime.dim[0] < 0 || hdr.dime.dim[0] > 15 ) )
		{
			hswf = true; // header swap
			hdr = headconv( hdr );
		}

		/**
		 * Assign field labeled with flag 2 (change) to current header.
		 */
		for ( unsigned int i = 0; i < niinfo; i++ )
		{
			if ( ( iinfo[i].flag & 2 ) == 2 )
			{
				assign( &hdr, iinfo[i] );
			}
		}

		/* re-swap header if necessary */
		if ( hswf )
		{
			hdr = headconv( hdr );
		}

		if ( zipped )
		{
			fdo = gzopen( inputFileName.c_str(), "w" );

			gzwrite( fdo, (char *) ( &hdr ), sizeof(struct dsr) );
			gzwrite( fdo, &( data[0] ), sizeof(unsigned char) * data.size() );
			gzclose( fdo );
		} else
		{
			if ( fseek( hdrfp, 0L, 0 ) != 0 )
			{
				std::cerr << "failed to rewind" << inputFileName << std::endl;
			}
			if ( fwrite( (char *) ( &hdr ), sizeof(struct dsr), 1, hdrfp ) != 1 )
				std::cerr << "failed to write" << inputFileName << std::endl;

			fclose( hdrfp );
		}
	}
};

/**
 * DBEDit.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::vector< std::string > inputs;
	float cal_min = std::numeric_limits< float >::min();
	float cal_max = std::numeric_limits< float >::max();
	int label = 0;

	tkd::CmdParser parser( argv[0], "Change image header properties." );

	parser.AddArgument( inputs, "inputs" )->AddAlias( "i" )->SetInput( "<strings>" )->SetDescription( "Input images" )->SetRequired( true )->SetMinMax(
			1, 10000 );

	parser.AddArgument( cal_min, "min" )->AddAlias( "min" )->SetInput( "<double>" )->SetDescription( "Minimal value" )->SetMinMax( 1, 1 );

	parser.AddArgument( cal_max, "max" )->AddAlias( "max" )->SetInput( "<double>" )->SetDescription( "Maximal value" )->SetMinMax( 1, 1 );

	parser.AddArgument( label, "label" )->AddAlias( "l" )->SetInput( "<int>" )->SetDescription( "Colormap: "
		"[0] -> 'Greyscale' "
		"[1] -> 'Red-Yellow' "
		"[2] -> 'Blue-Lightblue' "
		"[3] -> 'Red' "
		"[4] -> 'Green' "
		"[5] -> 'Blue' "
		"[6] -> 'Yellow' "
		"[7] -> 'Pink' "
		"[8] -> 'Cool' "
		"[9] -> 'Hot' "
		"[10] -> 'render1' "
		"[11] -> 'render2' "
		"[12] -> 'render3' "
		"[13] -> 'MGH Cortical' "
		"[14] -> 'MGH Subcortical' "
		"[15] -> 'Random-Rainbow' (default: 0)"

	)->SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	DBEdit dbedit = DBEdit();

	dbedit.Run( inputs, cal_min, cal_max, label );
}

