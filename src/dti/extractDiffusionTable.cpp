/*
 * Extract diffusion table and bvalue and write to text-file.
 *
 *  Created on: Jun 11, 2009
 *      Author: wim
 */

#include "tkdCmdParser.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"

#include "fidFID.h"
#include "fidFIDReader.h"
#include "fidFourier.h"

#include "tkdCmdParser.h"
#include "fidEPITools.h"
#include "fidCommon.h"

#include <iostream>
#include <fstream>

/**
 * Extract diffusion table and bvalues and write to output text files.
 */
class ExtractDiffusionTable
{

public:

	/**
	 * Extract
	 */
	void run( const std::string& fidFileName, const std::string& outputFileName )
	{

		fid::FIDReader::Pointer reader = fid::FIDReader::New();
		reader->SetFileName( fidFileName );

		try
		{
			reader->Read();
		} catch ( itk::ExceptionObject& e )
		{
			std::cerr << "Error reading reference FID: " << e.GetDescription() << std::endl;
			exit( EXIT_FAILURE );
		}

		fid::FID::Pointer fid = reader->GetFID();
		Procparser pp = fid->GetProcpar();

		int numberOfRO = pp.GetSize( "dro" );
		vnl_matrix< double > table( numberOfRO, 3 );

		for ( int i = 0; i < numberOfRO; ++i )
		{
			table( i, 0 ) = pp.GetAs< double > ( "dpe", i );
			table( i, 1 ) = pp.GetAs< double > ( "dro", i );
			table( i, 2 ) = pp.GetAs< double > ( "dsl", i );
		}

		/*
		 // both positive and negative
		 if ( pp.GetSize( "gdiff" ) == 1 ) {
		 vnl_matrix< double > twice( numberOfRO * 2, 3 );
		 for ( int i = 0; i < numberOfRO; ++i ) {
		 twice.set_row( i * 2, table.get_row( i ) );
		 twice.set_row( i * 2 + 1, table.get_row( i ) * -1. );
		 }
		 table = twice;
		 numberOfRO = 2 * numberOfRO;
		 }*/

		// Write to output.
		std::string orient = pp.GetAs< std::string > ( "orient" );

		if ( orient == "cor" )
		{
			vnl_matrix< double > permute( 3, 3 );
			permute.set_identity();

			permute( 0, 0 ) = -1;
			permute( 0, 1 ) = 0;
			permute( 0, 2 ) = 0;

			permute( 1, 0 ) = 0;
			permute( 1, 1 ) = 0;
			permute( 1, 2 ) = 1;

			permute( 2, 0 ) = 0;
			permute( 2, 1 ) = 1;
			permute( 2, 2 ) = 0;

			table = ( permute * table.transpose() ).transpose();
		} else if ( orient == "trans90" )
		{
			vnl_matrix< double > permute( 3, 3 );
			permute.set_identity();

			permute( 0, 0 ) = 0;
			permute( 0, 1 ) = 1;
			permute( 0, 2 ) = 0;

			permute( 1, 0 ) = -1;
			permute( 1, 1 ) = 0;
			permute( 1, 2 ) = 0;

			permute( 2, 0 ) = 0;
			permute( 2, 1 ) = 0;
			permute( 2, 2 ) = 1;

			table = ( permute * table.transpose() ).transpose();
		}

		//else if ( orient == "trans" || orient == "oblique") {
		else if ( orient == "trans" )
		{
			vnl_matrix< double > permute( 3, 3 );
			permute.set_identity();

			permute( 0, 0 ) = 1;
			permute( 0, 1 ) = 0;
			permute( 0, 2 ) = 0;

			permute( 1, 0 ) = 0;
			permute( 1, 1 ) = -1;
			permute( 1, 2 ) = 0;

			permute( 2, 0 ) = 0;
			permute( 2, 1 ) = 0;
			permute( 2, 2 ) = 1;

			table = ( permute * table.transpose() ).transpose();
		}

		else
		{
			std::cout << "WARNING: unsupported orientation: " << orient << std::endl;
		}

		double gdiff = pp.GetAs< double > ( "gdiff" );
		double tdelta = pp.GetAs< double > ( "tdelta" );
		double tDELTA = pp.GetAs< double > ( "tDELTA" );

		double bValue = ( 2. / vnl_math::pi ) * ( 2. / vnl_math::pi ) * ( tdelta * tdelta ) * ( gdiff * gdiff ) * ( tDELTA - tdelta / 3. )
				* ( 26751.98775 * 26751.98775 * 0.01 );

		// write output
		std::ofstream bvals( ( outputFileName + std::string( "_bvals" ) ).c_str() );
		std::ofstream bvecs( ( outputFileName + std::string( "_bvecs" ) ).c_str() );

		for ( int i = 0; i < numberOfRO; ++i )
		{
			vnl_vector< double > vector = table.get_row( i );
			bvals << ( vector( 0 ) * vector( 0 ) + vector( 1 ) * vector( 1 ) + vector( 2 ) * vector( 2 ) ) * bValue << std::endl;
			bvecs << vector( 0 ) << " " << vector( 1 ) << " " << vector( 2 ) << std::endl;
		}
	}
};

/**
 * Main.
 */
int main( int argc, char * argv[] )
{

	// arguments...
	std::string fid;
	std::string output;

	tkd::CmdParser parser( argv[0], "Extract diffusion table and b-values." );

	parser.AddArgument( fid, "fid" ) -> AddAlias( "f" ) -> SetInput( "<string>" ) -> SetDescription( "Varian Diffusion FID" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	parser.AddArgument( output, "outputroot" ) -> AddAlias( "o" ) -> SetInput( "<string>" ) -> SetDescription( "Output root" ) -> SetRequired(
			true ) -> SetMinMax( 1, 1 );

	if ( !parser.Parse( argc, argv ) )
	{
		parser.PrintUsage( std::cout );
		return EXIT_FAILURE;
	}

	ExtractDiffusionTable extractDiffusionTable = ExtractDiffusionTable();

	extractDiffusionTable.run( fid, output );
}

