namespace flens {

template <typename T>
void
copy(int N, const T *x, int incX, T *y, int incY)
{
    for (int i=0; i<N; ++i, x+=incX, y+=incY) {
        (*y) = (*x);
    }
}

template <typename T>
void
scal(const int N, T alpha, T *x, int incX)
{
    for (int i=0; i<N; ++i, x+=incX) {
        (*x) *= alpha;
    }
}

template <typename T>
void
axpy(int N, T alpha, const T *x, int incX, T *y, int incY)
{
    for (int i=0; i<N; ++i, x+=incX, y+=incY) {
        (*y) += alpha * (*x);
    }
}

template <typename T>
T
dot(int N, const T *x, int incX, const T *y, int incY)
{
    T result = T(0);
    for (int i=0; i<N; ++i, x+=incX, y+=incY) {
        result += (*x) * (*y);
    }
    return result;
}

} // namespace flens
