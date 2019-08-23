project( fftw3 )

include_directories( . )

add_library( fftwpp
Array.cc
Array.h
Complex.h
convolution.cc
convolution.h
fftw++.cc
fftw++.h
precision.h
utils.h )
