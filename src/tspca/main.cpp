#include "tspca.h"

std::string ConstructPath( const std::string& filename, const std::string& second, const std::string& extension )
{
  std::stringstream ss;
  ss << filename << second << "." << extension;
  return ss.str();
}

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "tspca", "Time series PCA on 4D image data" );

  std::string imageFileName;
  std::string maskFileName;
  std::string outputBase;
  int numberOfComponents = 0;
  float numberOfVariance = 0;

  p.AddArgument( imageFileName, "image" )
    ->AddAlias( "i" )
    ->SetDescription( "4D input image" )
    ->SetRequired( true );

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetDescription( "3D mask image" )
    ->SetRequired( true );

  p.AddArgument( outputBase, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Base output filename" )
    ->SetRequired( true );

  p.AddArgument( numberOfComponents, "components" )
    ->AddAlias( "nc" )
    ->SetDescription( "Retain number of components" );

  p.AddArgument( numberOfVariance, "variance" )
    ->AddAlias( "nv" )
    ->SetDescription( "Retain variance" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  PCA pca( numberOfComponents, numberOfVariance );
  pca.ReadImage( imageFileName );
  pca.ReadMask( maskFileName );
  pca.ImageToTimeSeries();
  pca.TransformPCA();

  pca.WriteDesignMatrix( ConstructPath( outputBase, "comp", "txt" ) );
  pca.WriteCoefficientMatrix( ConstructPath( outputBase, "coeff", "txt" ) );
  pca.WriteEigenvalues( ConstructPath( outputBase, "lambda", "txt" ) );
  pca.WriteReconstructedImage( ConstructPath( outputBase, "recon", "nii.gz" ) );
  pca.WriteRegressionCoefficients( ConstructPath( outputBase, "coeff", "nii.gz" ) );
  pca.WriteCorrelationImage( ConstructPath( outputBase, "cc", "nii.gz" ) );

  return 0;
}
