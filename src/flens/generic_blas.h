#ifndef FLENS_GENERIC_BLAS_H
#define FLENS_GENERIC_BLAS_H 1

namespace flens {

template <typename T>
    void
    copy(int N, const T *x, int incX, T *y, int incY);

template <typename T>
    void
    scal(const int N, T alpha, T *X, int incX);

template <typename T>
    void
    axpy(int N, T alpha, const T *y1, int incX, T *y2, int incY);

template <typename T>
    T
    dot(int N, const T *x, int incX, const T *Y, int incY);

} // namespace flens

#include <generic_blas.tcc>

#endif // FLENS_GENERIC_BLAS_H
