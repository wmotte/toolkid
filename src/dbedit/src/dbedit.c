#ifndef lint
static char sccsid[] = "@(#)dbedit++ 1.0 InvivoNMR (wim@invivonmr.uu.nl)";
#endif lint
#include <math.h>
#include "anautils.h"

#define mcalloc(a,b) calloc(a,b)

#define WORD -1
#define INT -2
#define FLOAT -3
#define CHAR -4
#define BYTE -5

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

#define hdrloc(_xyz_) ((int)(&(hdr._xyz_))-(int)(&hdr))

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
		"hist.smax", 0, INT, 0, "", 148 + 192, "hist.smin", 0, INT, 0, "", 148 + 196, };
#define niinfo sizeof(iinfo)/sizeof(Headinfo)

char *progname;
#define myerror(_a_, _b_) {char txt[1024]; (void)sprintf(txt, "%s: %s %s", progname, _a_, _b_); perror(txt); exit(1);}

/**
 * old main.
 */
int main_old( int argc, char** argv )
{
	int i, j, mode, hswf = 0;

	FILE *hdrfp;
	hdr.dime.dim[0] = 4;
	hdr.dime.dim[4] = 1;
	hdr.hk.regular = 'r';
	hdr.hk.sizeof_hdr = sizeof(struct dsr);
	progname = argv[0];

	if ( argc < 3 )
	{
		(void) fprintf( stderr, "usage: %s %s\n", argv[0], "file.hdr atr1[=val1]|all [atr2[=val2] [atr3[=val3] ..]]" );

		(void) fprintf( stderr, "attributes are:\n" );
		for ( i = 0; i < niinfo; i++ )
			(void) fprintf( stderr, "%s\n", iinfo[i].txt );
		(void) exit( 1 );
	}

	hdr.dime.dim[0] = 4;
	hdr.dime.dim[4] = 1;
	hdr.hk.regular = 'r';
	hdr.hk.sizeof_hdr = sizeof(struct dsr);

	for ( j = 2; j < argc; j++ )
	{
		char txt1[128], txt2[128];
		int nmatches;
		int match = -1;
		if ( j == 2 )
		{
			if ( !strcmp( argv[j], "all" ) )
			{
				for ( i = 0; i < niinfo; i++ )
					iinfo[i].flag |= 1;
				continue;
			}

		}
		if ( strlen( argv[j] ) > 127 )
		{
			(void) fprintf( stderr, "%s: %s - argument too long\n", argv[0], argv[j] );
			break;
		}
		nmatches = sscanf( argv[j], "%[].[a-zA-Z_0-9]=%s", txt1, txt2 );
		for ( i = 0; ( i < niinfo ) && ( match == -1 ); i++ )
			if ( !strcmp( txt1, iinfo[i].txt ) )
				match = i;
		if ( match == -1 )
		{
			(void) fprintf( stderr, "%s: cant decode %s - no matching attribute\n", argv[0], argv[j] );
			(void) exit( 1 );
		}

		if ( nmatches == 2 )
		{
			double atof();
			iinfo[match].flag |= 2;
			(void) strncpy( iinfo[match].valtxt, txt2, 63 );
			iinfo[match].val = (float) atof( iinfo[match].valtxt );
			mode = 1;
		} else if ( nmatches == 1 )
			iinfo[match].flag |= 1;
		else
		{
			(void) fprintf( stderr, "%s: cant decode %s\n", argv[0], argv[j] );
			exit( 1 );
		}
	}

	if ( mode == 1 )
	{
		if ( ( hdrfp = fopen( argv[1], "r+" ) ) == (FILE *) 0 )
		{
			if ( ( hdrfp = fopen( argv[1], "w" ) ) == (FILE *) 0 )
			{
				myerror("failed to open", argv[1]);
			}
		}
		else if ( fread( (char *) ( &hdr ), sizeof(struct dsr), 1, hdrfp ) != 1 )
		{
			(void) fprintf( stderr, "%s: warning: file %s too small\n", argv[0], argv[1] );
		}
	}
	else
	{
		if ( ( hdrfp = fopen( argv[1], "r" ) ) == (FILE *) 0 )
		{
			myerror("failed to open", argv[1]);
		}

		if ( fread( (char *) ( &hdr ), sizeof(struct dsr), 1, hdrfp ) != 1 )
		{
			(void) fprintf( stderr, "%s: warning: file %s too small\n", argv[0], argv[1] );
		}
	}

	/* check for swapped header */
	if ( ( hdr.hk.sizeof_hdr == 1543569408 ) || ( hdr.dime.dim[0] < 0 || hdr.dime.dim[0] > 15 ) )
	{
		hswf = 1;
		hdr = headconv( hdr );
	}

	for ( i = 0; i < niinfo; i++ )
	{
		if ( ( iinfo[i].flag & 2 ) == 2 )
		{
			assign( &hdr, iinfo[i] );
			mode = 1;
		}
		if ( ( iinfo[i].flag & 1 ) == 1 )
		{
			(void) printf( "%s\t", iinfo[i].txt );
			show( &hdr, iinfo[i] );
		}
	}

	switch ( hdr.dime.bitpix )
	{
	case 1:
		if ( hdr.dime.datatype != 1 )
			(void) fprintf( stderr,
			"%s: %s: warning: dime.datatype should be 1 (binary) when dime.bitpix is 1\n", argv[0], argv[1] );
		break;
	case 8: /* unsigned char */
		if ( hdr.dime.datatype != 2 )
			(void) fprintf( stderr,
			"%s: %s: warning: dime.datatype should be 2 (unsigned char) when dime.bitpix is 8\n", argv[0], argv[1] );
		break;
	case 16: /*    short      */
		if ( hdr.dime.datatype != 4 )
			(void) fprintf( stderr,
			"%s: %s: warning: dime.datatype should be 4 (short) when dime.bitpix is 16\n", argv[0], argv[1] );

		break;
	case 24: /*     rgb       */
		if ( hdr.dime.datatype != 128 )
			(void) fprintf( stderr,
			"%s: %s: warning: dime.datatype should be 128 (rgb) when dime.bitpix is 24\n", argv[0], argv[1] );
		break;
	case 32: /*    float  or int    */
		if ( hdr.dime.datatype != 16 && hdr.dime.datatype != 8 )
			(void) fprintf( stderr,
			"%s: %s: warning: dime.datatype should be 16 (float) or 8 (int) when dime.bitpix is 32\n", argv[0], argv[1] );
		break;
	case 64: /*   complex  or double   */
		if ( hdr.dime.datatype != 32 && hdr.dime.datatype != 64 )
			(void) fprintf( stderr,
			"%s: %s: warning: dime.datatype should be 32 (complex) or 64 (double) when dime.bitpix is 64\n", argv[0], argv[1] );
		break;
	default:
		(void) fprintf( stderr, "%s: %s: warning: dime.bitpix of %d unrecognised\n", argv[0], argv[1], hdr.dime.bitpix );
	}

	if ( mode == 1 )
	{
		if ( hswf )
		{ /* reswap header if necessary */
			hdr = headconv( hdr );
		}
		if ( fseek( hdrfp, 0L, 0 ) != 0 )myerror("failed to rewind", argv[1]);
		if ( fwrite( (char *) ( &hdr ), sizeof(struct dsr), 1, hdrfp ) != 1 )myerror("failed to write", argv[1]);
	}
	(void) fclose( hdrfp );
	return ( 0 );
}

assign( header, iinfo )
	struct dsr *header; Headinfo iinfo;
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
		(void) strncpy( &( hp[iinfo.loc] ), iinfo.valtxt, iinfo.type - 1 );
	}
}

show( header, iinfo )
	struct dsr *header;Headinfo iinfo;
{
	union
	{
		char b[4];
		short w;
		int i;
		float f;
	} un;
	char *hp = (char *) header;

	if ( iinfo.type == WORD )
	{
		un.b[0] = hp[iinfo.loc + 0];
		un.b[1] = hp[iinfo.loc + 1];
		(void) printf( "%d\n", (int) un.w );
	} else if ( iinfo.type == INT )
	{
		un.b[0] = hp[iinfo.loc + 0];
		un.b[1] = hp[iinfo.loc + 1];
		un.b[2] = hp[iinfo.loc + 2];
		un.b[3] = hp[iinfo.loc + 3];
		(void) printf( "%d\n", un.i );
	} else if ( iinfo.type == FLOAT )
	{
		un.b[0] = hp[iinfo.loc + 0];
		un.b[1] = hp[iinfo.loc + 1];
		un.b[2] = hp[iinfo.loc + 2];
		un.b[3] = hp[iinfo.loc + 3];
		(void) printf( "%g\n", un.f );
	} else if ( iinfo.type == CHAR )
		(void) printf( "%c\n", hp[iinfo.loc] );
	else if ( iinfo.type == BYTE )
		(void) printf( "%d\n", (int) hp[iinfo.loc] );
	else if ( iinfo.type > 0 )
		(void) printf( "\"%s\"\n", &( hp[iinfo.loc] ) );
}

