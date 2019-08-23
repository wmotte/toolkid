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

#ifndef FLENS_LISTINITIALIZER_H
#define FLENS_LISTINITIALIZER_H 1

namespace flens {

template <typename Container>
class ListInitializer
{
};

template <typename Engine>
    class DenseVector;

template <typename T>
    class Array;

template <typename T>
class ListInitializer<Array<T> >
{
    public:
        ListInitializer(T *begin, int length, T value);

        ListInitializer<Array<T> >
        operator,(T value);

    private:
        ListInitializer(T *begin, T *end, T value);

        T *_it;
        T *_end;
};

template <typename Engine>
class ListInitializer<DenseVector<Engine> >
{
    public:
        typedef typename Engine::ElementType T;

        ListInitializer(T *begin, int stride, int length, T value);

        ListInitializer<DenseVector<Engine> >
        operator,(T value);

        ListInitializer(T *begin, int stride, T *end, T value);

    private:
        T *_it;
        int _stride;
        T *_end;
};

template <typename Engine>
class GeMatrix;

template <typename Engine>
class ListInitializer<GeMatrix<Engine> >
{
    public:
        typedef typename Engine::ElementType T;

        ListInitializer(GeMatrix<Engine> &A, int row, int col, T value);

        ListInitializer<GeMatrix<Engine> >
        operator,(T value);

    private:
        GeMatrix<Engine> &_A;
        int _row, _col;
};

template <typename I>
class ListInitializerSwitch;

template <typename I>
class ListInitializerSwitch<DenseVector<I> >
{
    public:
        typedef typename I::ElementType T;

        ListInitializerSwitch(T *it, int stride, int length, T value);

        ListInitializerSwitch(const ListInitializerSwitch<DenseVector<I> > &l);

        ~ListInitializerSwitch();

        ListInitializer<DenseVector<I> >
        operator,(T x);

    private:
        T *_it;
        int _stride;
        T *_end;
        T _value;
        mutable bool _wipeOnDestruct;
};

template <typename I>
class ListInitializerSwitch<GeMatrix<I> >
{
    public:
        typedef typename I::ElementType T;

        ListInitializerSwitch(GeMatrix<I> &A, int row, int col, T value);

        ListInitializerSwitch(const ListInitializerSwitch<GeMatrix<I> > &l);

        ~ListInitializerSwitch();

        ListInitializer<GeMatrix<I> >
        operator,(T x);

    private:
        GeMatrix<I> &_A;
        int _row, _col;
        T _value;
        mutable bool _wipeOnDestruct;
};

} // namespace flens

#include <listinitializer.tcc>

#endif // FLENS_LISTINITIALIZER_H
