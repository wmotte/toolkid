#include "Complex.h"
#include "convolution.h"
#include "utils.h"

using namespace std;
using namespace fftwpp;

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// Number of iterations.
unsigned int N0=10000000;
unsigned int N=0;
unsigned int m=12;
unsigned int M=1;
  
bool Direct=false, Implicit=true, Explicit=false;

inline void init(Complex *e, Complex *f, Complex *g, unsigned int M=1) 
{
  unsigned int m1=m+1;
  unsigned int Mm=M*m1;
  double factor=1.0/cbrt((double) M);
  for(unsigned int i=0; i < Mm; i += m1) {
    double s=sqrt(1.0+i);
    double efactor=1.0/s*factor;
    double ffactor=(1.0+i)*s*factor;
    double gfactor=1.0/(1.0+i)*factor;
    Complex *ei=e+i;
    Complex *fi=f+i;
    Complex *gi=g+i;
    ei[0]=1.0*efactor;
    for(unsigned int k=1; k < m; k++) ei[k]=efactor*Complex(k,k+1);
    fi[0]=1.0*ffactor;
    for(unsigned int k=1; k < m; k++) fi[k]=ffactor*Complex(k,k+1);
    gi[0]=2.0*gfactor;
    for(unsigned int k=1; k < m; k++) gi[k]=gfactor*Complex(k,2*k+1);
  }
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
    int c = getopt(argc,argv,"hdeipM:N:m:n:");
    if (c == -1) break;
		
    switch (c) {
      case 0:
        break;
      case 'd':
        Direct=true;
        break;
      case 'e':
        Explicit=true;
        Implicit=false;
        break;
      case 'i':
        Implicit=true;
        Explicit=false;
        break;
      case 'p':
        break;
      case 'M':
        M=atoi(optarg);
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 'm':
        m=atoi(optarg);
        break;
      case 'n':
        N0=atoi(optarg);
        break;
      case 'h':
      default:
        usage(1);
    }
  }

  unsigned int n=4*m-3;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  if(N == 0) {
    N=N0/n;
    if(N < 10) N=10;
  }
  cout << "N=" << N << endl;
  
  unsigned int m1=m+1;
  unsigned int np=Explicit ? n/2+1 : m1;
  if(Implicit) np *= M;
    
  Complex *e=ComplexAlign(np);
  Complex *f=ComplexAlign(np);
  Complex *g=ComplexAlign(np);

  double *T=new double[N];

  if(Implicit) {
    ImplicitHBiConvolution C(m,M);
    Complex **E=new Complex *[M];
    Complex **F=new Complex *[M];
    Complex **G=new Complex *[M];
    for(unsigned int s=0; s < M; ++s) {
      unsigned int sm=s*m1;
      E[s]=e+sm;
      F[s]=f+sm;
      G[s]=g+sm;
    }
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g,M);
      seconds();
      C.convolve(E,F,G);
//      C.convolve(e,f,g);
      T[i]=seconds();
    }
    
    timings("Implicit",T,N);

    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << e[i] << endl;
    else cout << e[0] << endl;
  }
  
  if(Explicit) {
    ExplicitHBiConvolution C(n,m,f);
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g);
      seconds();
      C.convolve(e,f,g);
      T[i]=seconds();
    }
    
    timings("Explicit",T,N);

    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << e[i] << endl;
    else cout << e[0] << endl;
  }
  
  if(Direct) {
    DirectHBiConvolution C(m);
    init(e,f,g);
    Complex *h=ComplexAlign(m);
    seconds();
    C.convolve(h,e,f,g);
    T[0]=seconds();
    
    timings("Direct",T,1);

    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
    deleteAlign(h);
  }

  deleteAlign(g);
  deleteAlign(f);
  deleteAlign(e);
  
  delete [] T;
}
