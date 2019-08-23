#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "vnl_vector_to_std_vector.h"
#include "std_vector_to_vnl_vector.h"

#include <algorithm>
#include <iostream>

#include <vector>

namespace dki
{

	class Test
	{

	public:

		typedef double PixelType;

		typedef vnl_matrix< double > MatrixType;
		typedef vnl_vector< double > VectorType;

		typedef std::pair< double, vnl_vector< double > > EigenType;
		typedef std::vector< EigenType > EigenContainerType;

		/*
		 * Constructor.
		 */
		Test()
		{
			// **************************
			// Onces at init
			// **************************

			MatrixType index_table = GetIndexTable();
			MatrixType kvec_indices = GetKvecIndices();
			VectorType kvec_to_table = GetKvec2Table( GetRowSort( index_table ), kvec_indices );
			MatrixType indices_w0000 = GetIndicesW( 0, 0, 0, 0, index_table );
			MatrixType indices_w1111 = GetIndicesW( 1, 1, 1, 1, index_table );
			MatrixType indices_w2222 = GetIndicesW( 2, 2, 2, 2, index_table );
			MatrixType indices_w1122 = GetIndicesW( 1, 1, 2, 2, index_table );
			MatrixType indices_w0022 = GetIndicesW( 0, 0, 2, 2, index_table );
			MatrixType indices_w0011 = GetIndicesW( 0, 0, 1, 1, index_table );

			std::cout << "index_table: " << index_table << std::endl;
			std::cout << "kvec_indices: " << kvec_indices << std::endl;
			std::cout << "kvec_to_table: " << kvec_to_table << std::endl;

			// **************************
			// FOR each voxel
			// **************************

			// DT

			MatrixType DT( 3, 3 );

			DT( 0, 0 ) = 1;
			DT( 0, 1 ) = 0;
			DT( 0, 2 ) = 0.5;

			DT( 1, 0 ) = 0;
			DT( 1, 1 ) = 1;
			DT( 1, 2 ) = 0;

			DT( 2, 0 ) = 0.5;
			DT( 2, 1 ) = 0;
			DT( 2, 2 ) = 1;
			EigenContainerType eig = GetEigenSystem( DT );

			for ( unsigned int i = 0; i < eig.size(); i++ )
				std::cout << "eigenvalue: " << i << ": " << eig.at( i ).first << std::endl;

			for ( unsigned int i = 0; i < eig.size(); i++ )
				std::cout << "eigenvector: " << i << ": " << eig.at( i ).second << std::endl;

			// KT

			VectorType kvec( 15 );
			kvec( 0 ) = 1.1; kvec( 1 ) = 2.2; kvec( 2 ) = 3.3; kvec( 3 ) = 4.4;
			kvec( 4 ) = 5.5; kvec( 5 ) = 6.6; kvec( 6 ) = 7.7; kvec( 7 ) = 8.8;
			kvec( 8 ) = 9.9; kvec( 9 ) = 10.10; kvec( 10 ) = 11.11; kvec( 11 ) = 12.12;
			kvec( 12 ) = 13.13; kvec( 13 ) = 14.14; kvec( 14 ) = 15.15;

			VectorType kvec_extended = GetKvecExtended( kvec, kvec_to_table );

			std::cout << "kvec_extended: " << kvec_extended << std::endl;

			PixelType Wtilda0000 = GetWtilda( indices_w0000, eig, kvec_extended );
			PixelType Wtilda0011 = GetWtilda( indices_w0011, eig, kvec_extended );
			PixelType Wtilda0022 = GetWtilda( indices_w0022, eig, kvec_extended );
			PixelType Wtilda1111 = GetWtilda( indices_w1111, eig, kvec_extended );
			PixelType Wtilda1122 = GetWtilda( indices_w1122, eig, kvec_extended );
			PixelType Wtilda2222 = GetWtilda( indices_w2222, eig, kvec_extended );

			std::cout << "Wtilda1111: " << Wtilda0000 << std::endl;
			std::cout << "Wtilda1122: " << Wtilda0011 << std::endl;
			std::cout << "Wtilda1133: " << Wtilda0022 << std::endl;
			std::cout << "Wtilda2222: " << Wtilda1111 << std::endl;
			std::cout << "Wtilda2233: " << Wtilda1122 << std::endl;
			std::cout << "Wtilda3333: " << Wtilda2222 << std::endl;
		}

		/**
		 * Return eigenvalues and eigenvector sorted from l1 > l2 > l3.
		 */
		EigenContainerType GetEigenSystem( const MatrixType& DT )
		{
			vnl_symmetric_eigensystem< PixelType > eig( DT );

			EigenContainerType container( 3 );
			container.at( 0 ) = EigenType( eig.get_eigenvalue( 2 ), eig.get_eigenvector( 2 ) );
			container.at( 1 ) = EigenType( eig.get_eigenvalue( 1 ), eig.get_eigenvector( 1 ) );
			container.at( 2 ) = EigenType( eig.get_eigenvalue( 0 ), eig.get_eigenvector( 0 ) );
			return container;
		}
		/**
		 * Return sum eigenvalues.
		 */
		PixelType GetTrace( const MatrixType& DT )
		{
			return DT( 0, 0 ) + DT( 1, 1 ) + DT( 2, 2 );
		}

		/**
		 * Return fractional anisotropy (FA).
		 */
		PixelType GetFA( const MatrixType& DT )
		{
			PixelType isp = inner_product( DT, DT );
			if ( isp > 0.0 )
			{
				PixelType trace = GetTrace( DT );
				PixelType anisotropy = 3.0 * isp - trace * trace;
				PixelType fractionalAnisotropy = std::sqrt( anisotropy / ( 2.0 * isp ) );
				return fractionalAnisotropy;
			}
			return 0.0;
		}

		/**
		 * Index table W sort.
		 */
		MatrixType GetIndexTable()
		{
			MatrixType table( 81, 4 );

			unsigned int lin_idx = 0;

			for ( unsigned int i = 0; i < 3; i++ )
				for ( unsigned int j = 0; j < 3; j++ )
					for ( unsigned int k = 0; k < 3; k++ )
						for ( unsigned int l = 0; l < 3; l++ )
						{
							VectorType row( 4 );
							row( 0 ) = i;
							row( 1 ) = j;
							row( 2 ) = k;
							row( 3 ) = l;
							table.set_row( lin_idx, row );
							lin_idx++;
						}

			return table;
		}

		/**
		 * Sort rows in matrix ('ascend' way).
		 */
		MatrixType GetRowSort( const MatrixType& M )
		{
			MatrixType out( M );

			for ( unsigned int r = 0; r < out.rows(); r++ )
			{
				std::vector< PixelType > v = vnl_vector_to_std_vector( out.get_row( r ) );
				std::sort( v.begin(), v.end() );
				out.set_row( r, std_vector_to_vnl_vector( v ) );
			}

			return out;
		}

		/**
		 * Return matching indices.
		 */
		VectorType GetKvec2Table( const MatrixType& table, const MatrixType& kvec )
		{
			VectorType indices( table.rows(), 0 );

			for ( unsigned int i = 0; i < table.rows(); i++ )
			{
				for ( unsigned int j = 0; j < kvec.rows(); j++ )
				{
					if ( kvec.get_row( j ) == table.get_row( i ) )
					{
						indices( i ) = j;
					}
				}
			}
			return indices;
		}

		/**
		 * Kvec:
		 *
		 * k1111 k2222 k3333
		 * k1112 k1113 k1222
		 * k2223 k1333 k2333
		 * k1122 k1133 k2233
		 * k1123 k1223 k1233
		 */
		MatrixType GetKvecIndices()
		{
			MatrixType t( 15, 4 );

			PixelType data[60] =
			{ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 1, 1, 1, 1, 1, 2, 0, 2, 2, 2, 1, 2, 2, 2, 0, 0, 1, 1, 0, 0,
					2, 2, 1, 1, 2, 2, 0, 0, 1, 2, 0, 1, 1, 2, 0, 1, 2, 2 };

			t.set( data );

			return t;
		}

		/**
		 * Matlab sub2ind.
		 */
		unsigned int sub2ind( unsigned int nrow, unsigned int x, unsigned int y )
		{

			return ( x + y * nrow );
		}

		/**
		 * Return W_ijkl in terms of indices.
		 */
		MatrixType GetIndicesW( unsigned int i, unsigned int j, unsigned int k, unsigned int l, const MatrixType indices )
		{
			MatrixType W( 81, 4 );

			for ( unsigned int r = 0; r < W.rows(); r++ )
			{
				W( r, 0 ) = sub2ind( 3, indices( r, 0 ), i );
				W( r, 1 ) = sub2ind( 3, indices( r, 1 ), j );
				W( r, 2 ) = sub2ind( 3, indices( r, 2 ), k );
				W( r, 3 ) = sub2ind( 3, indices( r, 3 ), l );
			}

			return W;
		}

		/**
		 * Resort kvec.
		 */
		VectorType GetKvecExtended( const VectorType& kvec, const VectorType& indices )
		{
			VectorType out( indices.size(), 0 );

			for ( unsigned int i = 0; i < out.size(); i++ )
			{
				out( i ) = kvec( indices( i ) );
			}
			return out;
		}

		/**
		 * Get column-wise vector of sorted eigen vectors.
		 */
		VectorType EigenVectorsToStretchedFormat( const VectorType& e1, const VectorType& e2, const VectorType& e3 )
		{
			VectorType s( 9 );
			s( 0 ) = e1( 0 );
			s( 1 ) = e2( 0 );
			s( 2 ) = e3( 0 );

			s( 3 ) = e1( 1 );
			s( 4 ) = e2( 1 );
			s( 5 ) = e3( 1 );

			s( 6 ) = e1( 2 );
			s( 7 ) = e2( 2 );
			s( 8 ) = e3( 2 );

			return s;
		}

		/**
		 * Return W~, rotated in DT coordinates.
		 */
		PixelType GetWtilda( const MatrixType& indices, const EigenContainerType& eig, const VectorType& kvec )
		{
			VectorType V = EigenVectorsToStretchedFormat( eig.at( 0 ).second, eig.at( 1 ).second, eig.at( 2 ).second );

			VectorType col1( indices.rows() );
			VectorType col2( indices.rows() );
			VectorType col3( indices.rows() );
			VectorType col4( indices.rows() );

			for ( unsigned int i = 0; i < col1.size(); i++ )
			{
				col1( i ) = V( indices( i, 0 ) );
				col2( i ) = V( indices( i, 1 ) );
				col3( i ) = V( indices( i, 2 ) );
				col4( i ) = V( indices( i, 3 ) );
			}

			VectorType tmp = element_product( element_product( col1, col2 ), element_product( col3, col4 ) );

			return dot_product( tmp, kvec );
		}
	};
}

/**
 * Main
 */
int main()
{
	dki::Test();
	return 0;
}

