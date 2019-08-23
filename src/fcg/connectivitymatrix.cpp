#include "rsConnectivityMatrix.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "tkdCmdParser.h"

int main( int argc, char ** argv )
{
  tkd::CmdParser p( "connectivitymatrix", "Calculate ROI time series correlations" );

  std::string inputImageName;
  std::string labelImageName;
  std::string outputFileName;
  std::string outputLabelName;
  std::vector< int > skip;
  std::vector< int > labelRange;
  bool doPCA = false;

  p.AddArgument( inputImageName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input 4D time series image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output 2D matrix with r-values" )
    ->SetRequired( true );

  p.AddArgument( outputLabelName, "output-labels" )
    ->AddAlias( "ol" )
    ->SetDescription( "Output file for ROI labels" );

  p.AddArgument( labelImageName, "labels" )
    ->AddAlias( "l" )
    ->SetDescription( "Input 3D label image" )
    ->SetRequired( true );

  p.AddArgument( skip, "ignore" )
    ->AddAlias( "ig" )
    ->SetInput( "<integer integer>" )
    ->SetDescription( "Ignore number of time points at begin and end (default: 10 10)" )
    ->SetMinMax( 1, 2 );

  p.AddArgument( labelRange, "label-range" )
    ->AddAlias( "lr" )
    ->SetInput( "<integer integer>" )
    ->SetDescription( "Only include labels with values in [min, max]" )
    ->SetMinMax( 0, 2 );

  p.AddArgument( doPCA, "pca" )
    ->SetDescription( "Extract characteristic ROI profiles with PCA" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  typedef rs::ConnectivityMatrix::InputType InputType;
  typedef rs::ConnectivityMatrix::OutputType OutputType;
  typedef rs::ConnectivityMatrix::LabelType LabelType;
  typedef rs::ConnectivityMatrix::LabelSetType LabelSetType;

  typedef itk::ImageFileReader< InputType > InputReaderType;
  typedef itk::ImageFileReader< LabelType > LabelReaderType;
  typedef itk::ImageFileWriter< OutputType > WriterType;

  InputReaderType::Pointer inputReader = InputReaderType::New();
  inputReader->SetFileName( inputImageName.c_str() );
  inputReader->Update();

  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( labelImageName.c_str() );
  labelReader->Update();

  rs::ConnectivityMatrix cm;
  cm.SetInput( inputReader->GetOutput() );
  cm.SetLabels( labelReader->GetOutput() );
  cm.SetSkipTimePoints( skip.size() > 0 ? skip[ 0 ] : 10, skip.size() > 1 ? skip[ 1 ] : 10 );

  if ( labelRange.size() == 2 )
    {
    cm.SetLabelRange( labelRange[ 0 ], labelRange[ 1 ] );
    }

  cm.Run( doPCA );

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( cm.GetOutput() );
  writer->Update();

  if ( outputLabelName != "" )
    {
    std::ofstream out( outputLabelName.c_str() );
    LabelSetType labelSet = cm.GetLabelSet();
    for( LabelSetType::iterator i = labelSet.begin(); i != labelSet.end(); ++i )
      {
      out << i->first << std::endl;
      }
    }

  return 0;
}
