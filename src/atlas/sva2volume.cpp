/*
 * watershed.cpp
 *
 *  Created on: Jun 23, 2009
 *      Author: wim
 */
#include "itkImageRegionIteratorWithIndex.h"
#include "tkdCmdParser.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include <fstream>

#include "fidCommon.h"

typedef float PixelType;
typedef itk::Image< PixelType, 3 > ImageType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

/**
 * Convert Sparse Volume Format (.sva) Files to 3D volumes.
 *
 * These are plain-text files consisting of a two line header followed by multiple lines of commaseparated
 * coordinate data.
 * The header section consists of 2 colon-separated name/value pairs. The first is named
 * “Comment” the second is named “Dimensions”.
 *
 * For example, this header from a Smoothed Expression Energy volume file:
 *
 * Comment:Smoothed energy volume for gene Calb2 ImageseriesId 20619
 * Dimensions:67,41,58
 */

class Voxel {

private:
	PixelType m_value;
	ImageType::IndexType m_index;

public:
	Voxel( const ImageType::IndexType index, const PixelType value ) {
		m_index = index;
		m_value = value;
	}

	ImageType::IndexType getIndex() const {
		return m_index;
	}

	PixelType getValue() const {
		return m_value;
	}

	~Voxel() {
	}
};

class Sva2Volume {

public:

	/**
	 * Run.
	 */
	void run( const std::string& input, const std::string& output, double voxelSizeMm ) {

		createVolume( output, getOrigin( input, voxelSizeMm ), getSize( input ), getSpacing( voxelSizeMm ), getVoxels( input ) );
	}

protected:

	/**
	 * Create volume from given size and 4D pixelvectors. (x,y,z, id).
	 */
	void createVolume( const std::string& output, const ImageType::PointType& origin, const ImageType::SizeType& size,
			const ImageType::SpacingType& spacing, const std::vector< Voxel >& voxels ) {

		if ( voxels.empty() ) {
			std::cerr << "Voxels are empty" << std::endl;
			exit( EXIT_FAILURE );
		}

		ImageType::Pointer volume = ImageType::New();

		// set direction...
		ImageType::DirectionType direction;
		direction.SetIdentity();
		volume -> SetDirection( direction );

		// construct region...
		ImageType::RegionType region;
		region.SetSize( size );
		ImageType::IndexType index;
		index[0] = 0;
		index[1] = 0;
		index[2] = 0;

		region.SetIndex( index );
		volume -> SetOrigin( origin );
		volume -> SetSpacing( spacing );

		volume-> SetRegions( region );
		volume -> Allocate();

		for ( unsigned int i = 0; i < voxels.size(); i++ ) {
			Voxel voxel = voxels[i];
			ImageType::IndexType index = voxel.getIndex();
			PixelType value = voxel.getValue();
			volume -> SetPixel( index, value );
		}

		// write output...
		WriterType::Pointer writer = WriterType::New();

		// get volume in fslview position...
		ImageType::Pointer coronalVolume = reorientCoronal( volume );

		// get full volume ( flip half atlas -> full atlas )...
		ImageType::Pointer fullVolume = flipVolume( coronalVolume );

		writer -> SetInput( coronalVolume );
		writer -> SetFileName( output );

		try {
			writer -> Update();
		} catch ( itk::ExceptionObject& e ) {
			std::cerr << "Error writing: " << output << std::endl;
			std::cerr << e.GetDescription() << std::endl;
		}
	}

	/**
	 * Extract dims from sva file...
	 */
	ImageType::SizeType getSize( const std::string& input ) {
		ImageType::SizeType size;

		std::ifstream istr;
		istr.open( input.c_str() );
		if ( istr.fail() ) {
			std::cerr << "Failed to read: " << input << std::endl;
			exit( EXIT_FAILURE );
		} else {

			// get comment line...
			std::string comment;
			getline( istr, comment );

			// get dimensions...
			int dimx, dimy, dimz;
			std::string dimensions;
			getline( istr, dimensions );

			sscanf( dimensions.c_str(), "Dimensions:%d,%d,%d\n", &dimx, &dimy, &dimz );

			size[0] = dimx;
			size[1] = dimy;
			size[2] = dimz;
		}

		istr.close();
		return size;
	}

	/**
	 * Get Spacing.
	 */
	ImageType::SpacingType getSpacing( double mm ) {
		ImageType::SpacingType spacing;
		spacing[0] = mm;
		spacing[1] = mm;
		spacing[2] = mm;

		return spacing;
	}

	/**
	 * Extract origin from sva file, given input and voxel size in mm.
	 */
	ImageType::PointType getOrigin( const std::string& input, double mm ) {
		ImageType::PointType origin;

		std::ifstream istr;
		istr.open( input.c_str() );
		if ( istr.fail() ) {
			std::cerr << "Failed to read: " << input << std::endl;
			exit( EXIT_FAILURE );
		} else {

			// get comment line...
			std::string comment;
			getline( istr, comment );

			// get dimensions...
			int dimx, dimy, dimz;
			std::string dimensions;
			getline( istr, dimensions );

			sscanf( dimensions.c_str(), "Dimensions:%d,%d,%d\n", &dimx, &dimy, &dimz );

			origin[0] = ( dimx * mm * -0.5 );
			origin[1] = ( dimy * mm * -0.5 );
			origin[2] = ( dimz * mm * -0.5 );
		}

		istr.close();
		return origin;
	}

	/**
	 * get 2d vector, where first three vectors are x,y,z; last one is ID.
	 */
	std::vector< Voxel > getVoxels( const std::string& input ) {

		std::vector< Voxel > voxels;

		std::ifstream istr;
		istr.open( input.c_str() );
		if ( istr.fail() ) {
			std::cerr << "Failed to read: " << input << std::endl;
			exit( EXIT_FAILURE );
		} else {

			// skip first two lines...
			std::string remove;
			getline( istr, remove );
			getline( istr, remove );

			std::string line;

			while ( istr >> line ) {
				int x, y, z;
				PixelType value;
				sscanf( line.c_str(), "%d,%d,%d,%f\n", &x, &y, &z, &value );

				ImageType::IndexType index;
				index[0] = x;
				index[1] = y;
				index[2] = z;

				Voxel voxel( index, value );
				voxels.push_back( voxel );
			}
		}

		istr.close();

		return voxels;
	}

	/**
	 * Reset atlas to coronal slices.
	 */
	ImageType::Pointer reorientCoronal( ImageType::Pointer input ) {

		//bool flip[] = { false, true, true };
		bool flip[] = { true, true, true };
		int permutation[] = { 2, 1, 0 };

		input = fid::Common::Permute< PixelType, 3 >( fid::Common::Flip< PixelType, 3 >( input, flip, false ), permutation );

		ImageType::DirectionType direction;
		direction.SetIdentity();

		input -> SetDirection( direction );

		return input;
	}

	/**
	 * Flip all voxels != 0 to other size of volume.
	 */
	ImageType::Pointer flipVolume( ImageType::Pointer input ) {

		ImageType::SizeType inputSize = input -> GetLargestPossibleRegion().GetSize();

		unsigned int xsize = inputSize[0];

		IteratorType it( input, input -> GetLargestPossibleRegion() );

		for ( ; !it.IsAtEnd(); ++it ) {

			ImageType::IndexType index = it.GetIndex();
			PixelType intensity = input -> GetPixel( index );

			// if pixel > 0, flip to other side of volume center...
			if ( intensity > 0 ) {
				ImageType::IndexType newIndex = index;
				newIndex[0] = xsize - index[0];
				input -> SetPixel( newIndex, intensity );
			}
		}

		return input;
	}

};

/**
 * Main.
 */
int main( int argc, char * argv[] ) {

	// arguments...
	std::string input;
	std::string output;
	double voxelSizeMm = 0.1;

	tkd::CmdParser parser( argv[0], "Sparse Volume File (sva) To 3D Image Converter." );

	parser.AddArgument( input, "input" ) -> AddAlias( "i" ) -> SetInput( "<string>" ) -> SetDescription( "Input SVA File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "output" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output 3D File" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( voxelSizeMm, "voxel-size" ) -> AddAlias( "mm" ) -> SetInput( "<double>" ) -> SetDescription(
			"Voxel Size in mm (default: 0.1)" ) -> SetRequired( false ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) ) {
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	Sva2Volume sva2volume = Sva2Volume();

	sva2volume.run( input, output, voxelSizeMm );
}
