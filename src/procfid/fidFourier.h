#ifndef __fidFourier_h__
#define __fidFourier_h__
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_math.h"
#include "itkAutoPointer.h"

namespace fid
{

template< class TPrecision >
class Fourier
{
public:
	typedef TPrecision PrecisionType;
	typedef vcl_complex< PrecisionType > PixelType;
	typedef vnl_vector< PixelType > VectorType;
	typedef vnl_matrix< PixelType > MatrixType;

	Fourier( int N, int M = 0 );
	~Fourier();

	void SetOutputKSpace( bool );
	bool GetOutputKSpace() const;

	void ifft( VectorType& v );
	void ifft( MatrixType& m );

private:
	Fourier( const Fourier& );
	void operator=( const Fourier& );

	void ifft( VectorType& v, int n );

	struct Impl;
	itk::AutoPointer< Impl > m_Impl;
};

} // end namespace fid

#endif /*__fidFourier_h__*/
