#ifndef __ANAUTILS_H__
#define __ANAUTILS_H__
/* anautils.h
Header for anautils.c routines */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <float.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

#if defined(WIN32) && !defined(__CYGWIN32__)
#define FILESEP '\\'
#define EXTSEP '.'
#else
#define FILESEP '/'
#define EXTSEP '.'
#endif

/* Header define from Analyze authors */
#include "dbh.h"

#define A_HEADER_SIZE	348
#define BYTEBITS	8

struct dtypeinfo
{
  int numbits;
  int numels;
  int numbytes;
  double max;
  double min;
};

/* macros for data type byte swapping */
#define ICONV(a) (*(int *)cconv((char*) &(a), sizeof(int), 1))
#define FCONV(a) (*(float *)cconv((char*) &(a), sizeof(float), 1))
#define FCONVM(a,b) ((float *)cconv((char*) (a), sizeof(float), b))
#define SICONV(a) (*(short int *)cconv((char*) &(a), sizeof(short int), 1))
#define SICONVM(a,b) ((short int *)cconv((char*) (a), sizeof(short int), b))

/* function headers */
struct dtypeinfo get_datainfo( int type);
struct dsr headconv( struct dsr hdr );
char* ffext(char* path, char filesep, char extsep);
char *ffname(char *fname, const char* path, char filesep);
char *cconv(char *imgbuf, int elsize, size_t elnum);
int isdir(char * str);
int writehdr(char *iname, struct dsr hdr);
int writeimg(char *iname, char *img, int imgsize);
int writeana(char *iname, struct dsr hdr, char *img, int imgsize);

#endif
