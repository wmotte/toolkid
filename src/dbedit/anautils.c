/* anautils.c 
Utilities for manipulating analyze images
*/

#include "anautils.h"

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

  /* spm makes odd use of origin field */
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

struct dtypeinfo get_datainfo( int type)
{
    
    static int DataTypeTypes[9] = {0 ,1 ,2 ,4 ,8 ,16,32,64,128};
    static int DataTypeSizes[9] = {0 ,0 ,8 ,16,32,32,32,64,8  };
    static int DataTypeElnum[9] = {0 ,0 ,1 ,1 ,1 ,1 ,2 ,1 ,3  };
    static double DataTypeMin[9] = {0, 0, 0, -32768.0,-2147483648.0,
				    FLT_MIN,FLT_MIN,DBL_MIN,0  };
    static double DataTypeMax[9] = {0, 0, 255,32767, 2147483647.0,
				     FLT_MAX,FLT_MAX,DBL_MAX,255};

    int i;
    struct dtypeinfo dti;

    dti.numels=0;
    dti.numbits=0;
    dti.max=0;
    dti.min=0;
    for (i=1;i<=8;i++){
	if (type==DataTypeTypes[i]) {
	    dti.numbits = DataTypeSizes[i];
	    dti.numels = DataTypeElnum[i];
	    dti.max = DataTypeMax[i];
	    dti.min = DataTypeMin[i];
	    break;
	}
    }
    dti.numbytes=dti.numbits/BYTEBITS;
    return(dti);
}

char* ffext(char* path, char filesep, char extsep)
{
  char *i;
  char *fext = path+strlen(path);
  for (i=fext-1;i>=path;i--){
    if (*i==extsep){
      fext = i;
      break;
    }
    if (*i==filesep)
      break;
  }
  return(fext);
}

char* ffname(char* fname, const char* path, char filesep)
{
    char *tk=0, *ptk=0, *pcpy;
    char fileseps[2];

    sprintf(fileseps,"%c", filesep);
    pcpy = strdup(path);
    tk = strtok(pcpy,fileseps);
    while (tk!=0){
	ptk=tk;
	tk = strtok(0,fileseps);
    }
    strcpy(fname, ptk);
    free(pcpy);
    return(fname);
}

char* cconv(char *dbuf, int elsize, size_t elnum)
{
    int dd2, dm1;
    char *k,*m, *ep,tv;
    dd2 = (elsize+1) / 2;
    dm1 = elsize-1;
    k = dbuf;
    m = k+dm1;
    ep = dbuf+(elnum*elsize);
    while(k<ep){
	tv = *k;
	*k++ = *m;
	*m-- = tv;
	if (k>=m){
	    k+=dd2;
	    m=k+dm1;
	}	
    }
    return dbuf;
}

int isdir(char * str)
{
    struct stat stbuf;
    return ((stat(str, &stbuf) != -1) && ((stbuf.st_mode & S_IFMT) == S_IFDIR));
}

int writehdr(char *iname, struct dsr hdr)
{
  FILE * nfp;
  int i;
  strcpy(ffext(iname,FILESEP,EXTSEP), ".hdr");
  /* write header */
  nfp = fopen(iname, "wb");
  if (nfp == NULL) {
    return(0);
  }
  /* file ready to write */
  i = fwrite(&hdr, A_HEADER_SIZE, 1, nfp);
  fclose(nfp);
  return(1);
}

int writeimg(char *iname, char *img, int imgsize)
{
  FILE * nfp;
  int i;
  strcpy(ffext(iname,FILESEP,EXTSEP), ".img");
  /* write image */
  nfp = fopen(iname, "wb");
  if (nfp == NULL) {
    return(0);
  }
  i = fwrite(img,imgsize, 1, nfp);
  if (i!=1){
    fclose(nfp);
    return(0);
  }
  fclose(nfp);
  return(1);
}

int writeana(char *iname, struct dsr hdr, char *img, int imgsize)
{
  int i;
  if (i = writehdr(iname, hdr))
    i = writeimg(iname, img, imgsize);
  return(i);
}
