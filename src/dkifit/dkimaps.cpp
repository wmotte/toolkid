#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkDiffusionTensor3D.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>

#include "tkdCmdParser.h"

#include <algorithm> // for sort
#include <cmath> // for pow and sqrt
#include <valarray> // for atan
#include <iostream>
#include <iterator>
#include <vector>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "vnl_vector_to_std_vector.h"
#include "std_vector_to_vnl_vector.h"

namespace dki
{
	/**
	 * Calculate diffusion kurtosis maps from given input tensors.
	 *
	 * Modification to Tabesh et al: If eigenvalues are negative set them to 0!
	 *
	 */
	class DkiMaps
	{
	public:

		struct parameters
		{
			std::string inputFileName;
			std::string maskFileName;
			std::string outputFileName;
		};

		typedef double PixelType;

		typedef itk::Image< PixelType, 4 > ImageType;
		typedef itk::Image< PixelType, 3 > OutputImageType;
		typedef itk::ImageFileReader< OutputImageType > MaskReaderType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageLinearConstIteratorWithIndex< ImageType > ConstIterator4DType;
		typedef itk::ImageLinearIteratorWithIndex< ImageType > Iterator4DType;
		typedef itk::ImageRegionIteratorWithIndex< OutputImageType > Iterator3DType;

		typedef itk::ImageRegionConstIteratorWithIndex< OutputImageType > ConstIterator3DType;
		typedef itk::ImageFileWriter< OutputImageType > WriterType;

		typedef vnl_vector< PixelType > VectorType;
		typedef vnl_matrix< PixelType > MatrixType;

		typedef std::pair< PixelType, VectorType > EigenType;
		typedef std::vector< EigenType > EigenContainerType;

		/**
		 * Constructor.
		 */
		DkiMaps( const parameters& args )
		{
			SetInput( args.inputFileName );
			InitMask( args.maskFileName );
			AllocateOutput();
			InitGlobalIndices();
			CalculateMaps();
			Write( args.outputFileName );
		}

	protected:

		ImageType::Pointer m_Input;
		OutputImageType::Pointer m_Mask;


		// dti maps

		OutputImageType::Pointer m_FA;
		OutputImageType::Pointer m_Trace;
		OutputImageType::Pointer m_Lradial;
		OutputImageType::Pointer m_Laxial;

		// dki maps

		OutputImageType::Pointer m_MK;
		OutputImageType::Pointer m_Kradial;
		OutputImageType::Pointer m_Kaxial;

		// indices
		VectorType m_kvec_to_table;
		MatrixType m_indices_w0000;
		MatrixType m_indices_w1111;
		MatrixType m_indices_w2222;
		MatrixType m_indices_w1122;
		MatrixType m_indices_w0022;
		MatrixType m_indices_w0011;

		/**
		 * Init indices for Wtilda rotation.
		 */
		void InitGlobalIndices()
		{
			MatrixType index_table = GetIndexTable();
			MatrixType kvec_indices = GetKvecIndices();
			m_kvec_to_table = GetKvec2Table( GetRowSort( index_table ), kvec_indices );

			m_indices_w0000 = GetIndicesW( 0, 0, 0, 0, index_table );
			m_indices_w1111 = GetIndicesW( 1, 1, 1, 1, index_table );
			m_indices_w2222 = GetIndicesW( 2, 2, 2, 2, index_table );
			m_indices_w1122 = GetIndicesW( 1, 1, 2, 2, index_table );
			m_indices_w0022 = GetIndicesW( 0, 0, 2, 2, index_table );
			m_indices_w0011 = GetIndicesW( 0, 0, 1, 1, index_table );
		}

		/**
		 * DKI maps.
		 */
		void CalculateMaps()
		{
			ConstIterator4DType it( m_Input, m_Input->GetLargestPossibleRegion() );
			ConstIterator3DType mit( m_Mask, m_Mask->GetLargestPossibleRegion() );

			Iterator3DType itFA( m_FA, m_FA->GetLargestPossibleRegion() );
			Iterator3DType itTrace( m_Trace, m_Trace->GetLargestPossibleRegion() );
			Iterator3DType itLaxial( m_Laxial, m_Laxial->GetLargestPossibleRegion() );
			Iterator3DType itLradial( m_Lradial, m_Lradial->GetLargestPossibleRegion() );

			Iterator3DType itMK( m_MK, m_MK->GetLargestPossibleRegion() );
			Iterator3DType itKaxial( m_Kaxial, m_Kaxial->GetLargestPossibleRegion() );
			Iterator3DType itKradial( m_Kradial, m_Kradial->GetLargestPossibleRegion() );

			it.SetDirection( 3 );
			it.GoToBegin();

			mit.GoToBegin();

			itFA.GoToBegin();
			itTrace.GoToBegin();
			itLaxial.GoToBegin();
			itLradial.GoToBegin();

			itMK.GoToBegin();
			itKaxial.GoToBegin();
			itKradial.GoToBegin();

			if ( m_Input->GetLargestPossibleRegion().GetSize()[3] < 21 )
			{
				std::cerr << "*** ERROR ***: input does not contain less than 21 tensor elements!" << std::endl;
				exit( EXIT_FAILURE );
			}

			unsigned int sliceIndex = 0;

			// for each voxel

			while ( !it.IsAtEnd(), !mit.IsAtEnd() )
			{
				if ( mit.Get() != 0 )
				{

					// slice index
					if ( sliceIndex != it.GetIndex()[2] )
					{
						std::cout << "Processing slice: " << sliceIndex << std::endl;
						sliceIndex = it.GetIndex()[2];
					}

					VectorType b( 6 );
					VectorType kvec( 15 );

					while ( !it.IsAtEndOfLine() )
					{
						unsigned int i = ( it.GetIndex() )[3];

						if ( i < 6 )
							b( i ) = it.Get();
						else
							kvec( i - 6 ) = it.Get();

						++it;
					}

					MatrixType DT = GetDTITensor( b );

					itFA.Set( GetFA( DT ) );
					itTrace.Set( GetTrace( DT ) );

					EigenContainerType eig = GetEigenSystem( DT );

					// set negative eigenvalues to 0

					for ( unsigned int i = 0; i < eig.size(); i++ )
						if ( eig.at( i ).first < 0 )
							eig.at( i ).first = 0;

					itLaxial.Set( eig.at( 0 ).first );
					itLradial.Set( ( eig.at( 1 ).first + eig.at( 2 ).first ) / 2. );

					try
					{
						VectorType kvec_scaled = kvec * ( 1. / std::pow( GetTrace( DT ), 2 ) );
						VectorType k = GetKurtosisValues( kvec_scaled, eig ); // TODO
						itMK.Set( k( 0 ) );
						itKaxial.Set( k( 1 ) );
						itKradial.Set( k( 2 ) );
					} catch ( boost::math::evaluation_error e )
					{
						// convergence error, leave kurtosis output zero?
					}
				}

				// go to next voxel

				it.NextLine();
				++mit;
				++itFA;
				++itTrace;
				++itLaxial;
				++itLradial;

				++itMK;
				++itKaxial;
				++itKradial;
			}
		}

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
		PixelType GetRotatedW( const MatrixType& indices, const EigenContainerType& eig, const VectorType& kvec )
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

		/**
		 * Return MK, Kaxial and Kradial.
		 *
		 * Given eig
		 * Given kvec (15 kurtosis parameters)
		 */
		VectorType GetKurtosisValues( const VectorType kvec, const EigenContainerType& eig )
		{
			VectorType kvec_extended = GetKvecExtended( kvec, m_kvec_to_table );

			PixelType Wtilda_0000 = GetRotatedW( m_indices_w0000, eig, kvec_extended );
			PixelType Wtilda_0011 = GetRotatedW( m_indices_w0011, eig, kvec_extended );
			PixelType Wtilda_0022 = GetRotatedW( m_indices_w0022, eig, kvec_extended );
			PixelType Wtilda_1111 = GetRotatedW( m_indices_w1111, eig, kvec_extended );
			PixelType Wtilda_1122 = GetRotatedW( m_indices_w1122, eig, kvec_extended );
			PixelType Wtilda_2222 = GetRotatedW( m_indices_w2222, eig, kvec_extended );

			PixelType l1 = eig.at( 0 ).first;
			PixelType l2 = eig.at( 1 ).first;
			PixelType l3 = eig.at( 2 ).first;

			PixelType MK = F1( l1, l2, l3 ) * Wtilda_0000 + F1( l2, l1, l3 ) * Wtilda_1111 + F1( l3, l2, l1 ) * Wtilda_2222 + F2( l1, l2,
					l3 ) * Wtilda_1122 + F2( l2, l1, l3 ) * Wtilda_0022 + F2( l3, l2, l1 ) * Wtilda_0011;

			PixelType Kaxial = ( std::pow( l1 + l2 + l3, 2 ) / 9 * std::pow( l1, 2 ) ) * Wtilda_0000;

			PixelType Kradial = G1( l1, l2, l3 ) * Wtilda_1111 + G1( l1, l3, l2 ) * Wtilda_2222 + G2( l1, l2, l3 ) * Wtilda_1122;

			VectorType k( 3 );
			k( 0 ) = MK;
			k( 1 ) = Kaxial;
			k( 2 ) = Kradial;

			return k;
		}

		/**
		 * Eq (33) Tabesh et al.
		 */
		PixelType G1( PixelType l1, PixelType l2, PixelType l3 )
		{
			if ( l2 != l3 )
			{
				PixelType a = std::pow( l1 + l2 + l3, 2 ) / 18 * l1 * std::pow( l2 - l3, 2 );
				PixelType b = 2 * l2 + ( std::pow( l3, 2 ) - 3 * l2 * l3 ) / std::sqrt( l2 * l3 );

				return a * ( 2 * l2 + b );
			}

			return std::pow( l1 + 2 * l2, 2 ) / std::pow( 24. * l2, 2 );
		}

		/**
		 * Eq (34) Tabesh et al.
		 */
		PixelType G2( PixelType l1, PixelType l2, PixelType l3 )
		{
			if ( l2 != l3 )
			{
				PixelType a = std::pow( l1 + l2 + l3, 2 ) / 3 * std::pow( l2 - l3, 2 );
				PixelType b = l2 + l3 / std::sqrt( l2 + l3 );

				return a * ( b - 2 );
			}

			return 6 * std::pow( l1 + 2 * l2, 2 ) / std::pow( 72 * l2, 2 );
		}

		/**
		 * Eq (27) Tabesh et al.
		 */
		PixelType F1( PixelType l1, PixelType l2, PixelType l3 )
		{
			if ( ( l1 == l2 ) && ( l2 == l3 ) )
			{
				return 0.2;
			} else if ( l1 == l2 )
			{
				return 0.5 * F2( l3, l1, l1 );
			} else if ( l1 == l3 )
			{
				return 0.5 * F2( l2, l1, l1 );
			}

			PixelType a = std::pow( l1 + l2 + l3, 2 ) / 18 * ( l1 - l2 ) * ( l1 - l3 );
			PixelType b = std::sqrt( l2 * l3 ) / l1;
			PixelType c = boost::math::ellint_rf( l1 / l2, l1 / l3, 1. );
			PixelType d = ( 3 * std::pow( l1, 2 ) - l1 * l2 - l1 * l3 - l2 * l3 ) / 3 * l1 * std::sqrt( l2 * l3 );
			PixelType e = boost::math::ellint_rd( l1 / l2, l1 / l3, 1. );

			return a * ( b * c + d * e - 1 );
		}

		/**
		 * Eq (28) Tabesh et al.
		 */
		PixelType F2( PixelType l1, PixelType l2, PixelType l3 )
		{
			if ( l2 != l3 )
			{
				PixelType a = std::pow( l1 + l2 + l3, 2 ) / ( 3 * std::pow( l2 - l3, 2 ) );
				PixelType b = ( l2 + l3 ) / std::sqrt( l2 * l3 );
				PixelType c = boost::math::ellint_rf( l1 / l2, l1 / l3, 1. );
				PixelType d = ( 2 * l1 - l2 + l3 ) / ( 3 * std::sqrt( l2 * l3 ) );
				PixelType e = boost::math::ellint_rd( l1 / l2, l1 / l3, 1. );

				return a * ( b * c + d * e - 2 );
			} else if ( l1 != l3 )
			{
				PixelType a = std::pow( l1 + 2 * l3, 2 ) / std::pow( 144 * l3, 2 ) * std::pow( l1 - l3, 2 );
				PixelType b = l3 * ( l1 + 2 * l3 );
				PixelType c = 1 - ( l1 / l3 );
				PixelType d = l1 * ( l1 - 4 * l3 );

				c = std::sqrt( std::abs< PixelType >( c ) );
				d = d * ( 1. / c ) * std::atan( c );

				return 6 * a * ( b + d );
			}

			return 6. / 15.;
		}

		/**
		 * Return 2-rank DTI tensor.
		 *
		 * layout =>
		 * 		| 0  1  2  |
		 *	    | 1  3  4  |
		 *      | 2  4  5  |
		 */
		MatrixType GetDTITensor( const VectorType& x )
		{
			MatrixType DT( 3, 3 );

			/*
			 m_A_D( i, 0 ) =     Gx * Gx;
			 m_A_D( i, 1 ) =     Gy * Gy;
			 m_A_D( i, 2 ) = 2 * Gx * Gy;
			 m_A_D( i, 3 ) =     Gz * Gz;
			 m_A_D( i, 4 ) = 2 * Gy * Gz;
			 m_A_D( i, 5 ) = 2 * Gx * Gz;
			 */

			DT( 0, 0 ) = x( 0 ); // x * x
			DT( 0, 1 ) = x( 2 );
			DT( 0, 2 ) = x( 5 );

			DT( 1, 0 ) = x( 2 );
			DT( 1, 1 ) = x( 1 ); // y * y
			DT( 1, 2 ) = x( 4 );

			DT( 2, 0 ) = x( 5 );
			DT( 2, 1 ) = x( 4 );
			DT( 2, 2 ) = x( 3 ); // z * z

			return DT;
		}

		/**
		 * Set input image.
		 */
		void SetInput( const std::string& inputFileName )
		{
			if ( !inputFileName.empty() )
			{
				ReaderType::Pointer reader = ReaderType::New();
				reader->SetFileName( inputFileName );
				reader->Update();
				m_Input = reader->GetOutput();
			} else
			{
				std::cerr << "Could not read input: " << inputFileName << "!" << std::endl;
				exit( EXIT_FAILURE );
			}
		}

		/**
			 * Init mask if file given, else create empty mask from input file.
			 */
			void InitMask( const std::string& maskFileName )
			{
				if ( !maskFileName.empty() )
				{
					MaskReaderType::Pointer reader = MaskReaderType::New();
					reader->SetFileName( maskFileName );
					reader->Update();
					m_Mask = reader->GetOutput();
				} else
				{
					OutputImageType::Pointer output = OutputImageType::New();
					OutputImageType::RegionType region3D;
					OutputImageType::IndexType index3D;
					OutputImageType::SizeType size3D;
					OutputImageType::SpacingType spacing3D;
					OutputImageType::PointType origin3D;

					ImageType::RegionType region = m_Input->GetLargestPossibleRegion();
					ImageType::SizeType size = region.GetSize();
					ImageType::IndexType index = region.GetIndex();
					ImageType::SpacingType spacing = m_Input->GetSpacing();
					ImageType::PointType origin = m_Input->GetOrigin();

					size3D[0] = size[0];
					size3D[1] = size[1];
					size3D[2] = size[2];
					index3D[0] = index[0];
					index3D[1] = index[1];
					index3D[2] = index[2];
					origin3D[0] = origin[0];
					origin3D[1] = origin[1];
					origin3D[2] = origin[2];
					spacing3D[0] = spacing[0];
					spacing3D[1] = spacing[1];
					spacing3D[2] = spacing[2];

					region3D.SetSize( size3D );
					region3D.SetIndex( index3D );

					// set
					output->SetRegions( region3D );
					output->SetSpacing( spacing3D );
					output->SetOrigin( origin3D );
					output->Allocate();
					output->FillBuffer( 1 );

					m_Mask = output;
				}
			}

		/**
		 * Allocate all output images to 0.
		 */
		void AllocateOutput()
		{
			OutputImageType::RegionType region3D;
			OutputImageType::IndexType index3D;
			OutputImageType::SizeType size3D;
			OutputImageType::SpacingType spacing3D;
			OutputImageType::PointType origin3D;

			ImageType::RegionType region = m_Input->GetLargestPossibleRegion();
			ImageType::SizeType size = region.GetSize();
			ImageType::IndexType index = region.GetIndex();
			ImageType::SpacingType spacing = m_Input->GetSpacing();
			ImageType::PointType origin = m_Input->GetOrigin();

			size3D[0] = size[0];
			size3D[1] = size[1];
			size3D[2] = size[2];
			index3D[0] = index[0];
			index3D[1] = index[1];
			index3D[2] = index[2];
			origin3D[0] = origin[0];
			origin3D[1] = origin[1];
			origin3D[2] = origin[2];
			spacing3D[0] = spacing[0];
			spacing3D[1] = spacing[1];
			spacing3D[2] = spacing[2];

			region3D.SetSize( size3D );
			region3D.SetIndex( index3D );

			// create
			m_FA = OutputImageType::New();
			m_Trace = OutputImageType::New();
			m_Laxial = OutputImageType::New();
			m_Lradial = OutputImageType::New();

			m_MK = OutputImageType::New();
			m_Kaxial = OutputImageType::New();
			m_Kradial = OutputImageType::New();

			// set
			m_FA->SetRegions( region3D );
			m_FA->SetSpacing( spacing3D );
			m_FA->SetOrigin( origin3D );
			m_FA->Allocate();
			m_FA->FillBuffer( 0 );

			m_Trace->SetRegions( region3D );
			m_Trace->SetSpacing( spacing3D );
			m_Trace->SetOrigin( origin3D );
			m_Trace->Allocate();
			m_Trace->FillBuffer( 0 );

			m_Laxial->SetRegions( region3D );
			m_Laxial->SetSpacing( spacing3D );
			m_Laxial->SetOrigin( origin3D );
			m_Laxial->Allocate();
			m_Laxial->FillBuffer( 0 );

			m_Lradial->SetRegions( region3D );
			m_Lradial->SetSpacing( spacing3D );
			m_Lradial->SetOrigin( origin3D );
			m_Lradial->Allocate();
			m_Lradial->FillBuffer( 0 );

			m_MK->SetRegions( region3D );
			m_MK->SetSpacing( spacing3D );
			m_MK->SetOrigin( origin3D );
			m_MK->Allocate();
			m_MK->FillBuffer( 0 );

			m_Kaxial->SetRegions( region3D );
			m_Kaxial->SetSpacing( spacing3D );
			m_Kaxial->SetOrigin( origin3D );
			m_Kaxial->Allocate();
			m_Kaxial->FillBuffer( 0 );

			m_Kradial->SetRegions( region3D );
			m_Kradial->SetSpacing( spacing3D );
			m_Kradial->SetOrigin( origin3D );
			m_Kradial->Allocate();
			m_Kradial->FillBuffer( 0 );
		}

		/**
		 * Write output images.
		 */
		void Write( const std::string& outputFileName )
		{
			// Writer
			WriterType::Pointer writer = WriterType::New();

			// FA
			writer->SetFileName( ( outputFileName + std::string( "_FA.nii.gz" ) ).c_str() );
			writer->SetInput( m_FA );
			writer->Update();

			// Trace
			writer->SetFileName( ( outputFileName + std::string( "_Trace.nii.gz" ) ).c_str() );
			writer->SetInput( m_Trace );
			writer->Update();

			// L axial
			writer->SetFileName( ( outputFileName + std::string( "_Laxial.nii.gz" ) ).c_str() );
			writer->SetInput( m_Laxial );
			writer->Update();

			// L radial
			writer->SetFileName( ( outputFileName + std::string( "_Lradial.nii.gz" ) ).c_str() );
			writer->SetInput( m_Lradial );
			writer->Update();

			// MK
			writer->SetFileName( ( outputFileName + std::string( "_MK.nii.gz" ) ).c_str() );
			writer->SetInput( m_MK );
			writer->Update();

			// K axial
			writer->SetFileName( ( outputFileName + std::string( "_Kaxial.nii.gz" ) ).c_str() );
			writer->SetInput( m_Kaxial );
			writer->Update();

			// K radial
			writer->SetFileName( ( outputFileName + std::string( "_Kradial.nii.gz" ) ).c_str() );
			writer->SetInput( m_Kradial );
			writer->Update();
		}

	};

} // end namespace dki


/**
 * Create dki maps.
 */
int main( int argc, char ** argv )
{
	tkd::CmdParser p( argv[0], "Create diffusion kurtosis maps." );

	dki::DkiMaps::parameters args;

	p.AddArgument( args.inputFileName, "input" ) ->AddAlias( "i" ) ->SetDescription( "Input 4D image with 21 tensor elements" ) ->SetRequired(
			true );

	p.AddArgument( args.outputFileName, "output" ) ->AddAlias( "o" ) ->SetDescription( "Output filename base" ) ->SetRequired( true );

	p.AddArgument( args.maskFileName, "mask" ) ->AddAlias( "m" ) ->SetDescription( "Mask 3D image" );

	if ( !p.Parse( argc, argv ) )
	{
		p.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	dki::DkiMaps maps( args );

	return EXIT_SUCCESS;
}

