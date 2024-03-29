#include "Complex.h"
#include "convolution.h"
#include "utils.h"
#include "Array.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv2.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv2.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// Number of iterations.
unsigned int N0=10000000;
unsigned int N=0;
unsigned int nx=0;
unsigned int ny=0;
unsigned int mx=4;
unsigned int my=4;
unsigned int M=1;

unsigned int nxp;
unsigned int nyp;

bool Direct=false, Implicit=true, Explicit=false, Pruned=false;

unsigned int outlimit=100;

inline void init(array2<Complex>& e, array2<Complex>& f, array2<Complex>& g,
                 unsigned int M=1) 
{
  unsigned int offset=Explicit ? nx/2-mx+1 : (Implicit ? 1 : 0);
  unsigned int stop=2*mx-1;
  unsigned int stopoffset=stop+offset;
  double factor=1.0/cbrt((double) M);
  for(unsigned int s=0; s < M; ++s) {
    double S=sqrt(1.0+s);
    double efactor=1.0/S*factor;
    double ffactor=(1.0+S)*S*factor;
    double gfactor=1.0/(1.0+S)*factor;
    for(unsigned int i=0; i < stop; i++) {
      for(unsigned int j=0; j < my; j++) {
        unsigned int I=s*stopoffset+i+offset;
        e[I][j]=efactor*Complex(i,j);
        f[I][j]=ffactor*Complex(i+1,j+2);
        g[I][j]=gfactor*Complex(2*i,j+1);
      }
    }
  }
}

unsigned int padding(unsigned int m)
{
  unsigned int n=4*m-3;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  return 1 << log2n;
}

int main(int argc, char* argv[])
{
#ifndef __SSE2__
  fftw::effort |= FFTW_NO_SIMD;
#endif  
  
#ifdef __GNUC__	
  optind=0;
#endif	
  for (;;) {
    int c = getopt(argc,argv,"hdeiptM:N:m:x:y:n:");
    if (c == -1) break;
		
    switch (c) {
      case 0:
        break;
      case 'd':
        Direct=true;
        Implicit=false;
        break;
      case 'e':
        Explicit=true;
        Implicit=false;
        Pruned=false;
        break;
      case 'i':
        Implicit=true;
        Direct=false;
        Explicit=false;
        break;
      case 'p':
        Explicit=true;
        Implicit=false;
        Pruned=true;
        break;
      case 'M':
        M=atoi(optarg);
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 'm':
        mx=my=atoi(optarg);
        break;
      case 'x':
        mx=atoi(optarg);
        break;
      case 'y':
        my=atoi(optarg);
        break;
      case 'n':
        N0=atoi(optarg);
        break;
      case 'h':
      default:
        usage(2);
    }
  }

  nx=padding(mx);
  ny=padding(my);
  
  cout << "nx=" << nx << ", ny=" << ny << endl;
  cout << "mx=" << mx << ", my=" << my << endl;
  
  if(N == 0) {
    N=N0/(nx*ny);
    if(N < 10) N=10;
  }
  cout << "N=" << N << endl;
    
  size_t align=sizeof(Complex);
  nxp=Explicit ? nx : (Implicit ? 2*mx : 2*mx-1);
  nyp=Explicit ? ny/2+1 : (Implicit ? my+1 : my);
  unsigned int nxp0=Implicit ? nxp*M : nxp;
  array2<Complex> e(nxp0,nyp,align);
  array2<Complex> f(nxp0,nyp,align);
  array2<Complex> g(nxp0,nyp,align);

  double *T=new double[N];

  if(Implicit) {
    ImplicitHBiConvolution2 C(mx,my,M);
    Complex **E=new Complex *[M];
    Complex **F=new Complex *[M];
    Complex **G=new Complex *[M];
    unsigned int mf=nxp*nyp;
    for(unsigned int s=0; s < M; ++s) {
      unsigned int smf=s*mf;
      E[s]=e+smf;
      F[s]=f+smf;
      G[s]=g+smf;
    }
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g,M);
      seconds();
      C.convolve(E,F,G);
//      C.convolve(e,f,g);
      T[i]=seconds();
    }
    
    timings("Implicit",T,N);
    
    if(nxp*my < outlimit)
      for(unsigned int i=1; i < nxp; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << e[i][j] << "\t";
        cout << endl;
      } else cout << e[1][0] << endl;
    cout << endl;
  }
  
  if(Explicit) {
    ExplicitHBiConvolution2 C(nx,ny,mx,my,f,Pruned);
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g);
      seconds();
      C.convolve(e,f,g);
      T[i]=seconds();
    }
    
    timings(Pruned ? "Pruned" : "Explicit",T,N);

    unsigned int offset=nx/2-mx+1;
    if(2*(mx-1)*my < outlimit) 
      for(unsigned int i=offset; i < offset+2*mx-1; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << e[i][j] << "\t";
        cout << endl;
      } else cout << e[offset][0] << endl;
  }
  
  if(Direct) {
    Explicit=0;
    unsigned int nxp=2*mx-1;
    array2<Complex> h(nxp,my,align);
    array2<Complex> e(nxp,my,align);
    array2<Complex> f(nxp,my,align);
    array2<Complex> g(nxp,my,align);
    DirectHBiConvolution2 C(mx,my);
    init(e,f,g);
    seconds();
    C.convolve(h,e,f,g);
    T[0]=seconds();
  
    timings("Direct",T,1);

    if(nxp*my < outlimit)
      for(unsigned int i=0; i < nxp; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << h[i][j] << "\t";
        cout << endl;
      } else cout << h[0][0] << endl;
  }
  
  delete [] T;
}
