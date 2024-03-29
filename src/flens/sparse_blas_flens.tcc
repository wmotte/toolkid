/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cassert>
#include <sparse_blas.h>
#include <sparse_blas_flens.h>

namespace flens {

template <typename T, typename VX, typename VY>
void
mv(Transpose trans, T alpha, const SparseGeMatrix<CRS<T> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::T beta, DenseVector<VY> &y)
{
    assert(ADDRESS(y)!=ADDRESS(x));
    if (trans==NoTrans) {
        assert(A.numCols()==x.length());
    }
    if (trans==Trans) {
        assert(A.numRows()==x.length());
    }
    int yLength = (trans==NoTrans) ? A.numRows()
                                   : A.numCols();

    assert((beta==0) || (y.length()==yLength));

    if (y.length()!=yLength) {
        y.resize(yLength);
    }

    assert(x.stride()==1);
    assert(y.stride()==1);

    crs_gemv(trans, A.numRows(), A.numCols(),
             alpha,
             A.engine().values.data(),
             A.engine().rows.data(),
             A.engine().columns.data(),
             x.data(),
             beta, y.data());
}

// sparse_symv
template <typename T, CRS_Storage Storage, typename VX, typename VY>
void
mv(T alpha, const SparseSyMatrix<CRS<T, Storage> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::T beta, DenseVector<VY> &y)
{
    assert(A.dim()==x.length());
    assert(A.dim()==y.length());

    assert(x.stride()==1);
    assert(y.stride()==1);

    if (Storage==CRS_General) {
        crs_gemv(NoTrans, A.dim(), A.dim(),
                 alpha,
                 A.engine().values.data(),
                 A.engine().rows.data(),
                 A.engine().columns.data(),
                 x.data(),
                 beta, y.data());
    } else {
        StorageUpLo upLo = (Storage==CRS_UpperTriangular) ? Upper : Lower;
        crs_symv(upLo, A.dim(), alpha,
                 A.engine().values.data(),
                 A.engine().rows.data(),
                 A.engine().columns.data(),
                 x.data(),
                 beta, y.data());
    }
}

// sparse_gemm (matrix B dense!)
template <typename T, CRS_Storage Storage>
void
mm(Transpose transA,
   Transpose transB,
   T alpha,
   const SparseGeMatrix<CRS<T, Storage> > &A,
   const GeMatrix<FullStorage<T, ColMajor> > &B,
   T beta, GeMatrix<FullStorage<T, ColMajor> > &C)
{
    assert(Storage==CRS_General);
    assert(transB==NoTrans); // Transposition not provided for B.
    int m = (transA==NoTrans) ? A.numRows() : A.numCols();
    int n = B.numCols();

    assert((beta==0) || (C.numRows()==m));
    assert((beta==0) || (C.numCols()==n));

    C.resize(m,n);

    int r = A.numRows();
    DenseVector<Array<int> > pointerE(r,1);
    pointerE(_(1,r-1)) = A.engine().rows(_(2,r));
    pointerE(r) = A.engine().values.length()+1;

    char matdescra[] = {'G', 'U', 'N', 'F', ' ', ' '};

    csrmm(transA, A.numRows(), C.numCols(), A.numCols(), alpha, matdescra,
          A.engine().values.data(), A.engine().columns.data(),
          A.engine().rows.data(), pointerE.data(),
          B.data(), B.leadingDimension(), beta, C.data(), C.leadingDimension());
}

} // namespase flens
