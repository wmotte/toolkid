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

#ifndef FLENS_FLENS_H
#define FLENS_FLENS_H 1

#define ADDRESS(x) reinterpret_cast<const void *>(&x)

#ifndef ASSERT
#define ASSERT(x) assert(x)
#endif //ASSERT

#include <array.h>
#include <bandstorage.h>
#include <blas.h>
#include <blas_flens.h>
#include <cg.h>
#include <complex_helper.h>
#include <crs.h>
#include <densevector.h>
#include <evalclosure.h>
#include <fullstorage.h>
#include <generalmatrix.h>
#include <generic_blas.h>
#include <refcounter.h>
#include <lapack.h>
#include <lapack_flens.h>
#include <listinitializer.h>
#include <matvec.h>
#include <matvecclosures.h>
#include <matvecio.h>
#include <matvecoperations.h>
#include <multigrid.h>
#include <packedstorage.h>
#include <range.h>
#include <snapshot.h>
#include <storage.h>
#include <sparsematrix.h>
#include <sparse_blas.h>
#include <sparse_blas_flens.h>
#include <symmetricmatrix.h>
#include <traits.h>
#include <triangularmatrix.h>
#include <underscore.h>
#include <uplo.h>

#endif // FLENS_FLENS_H
