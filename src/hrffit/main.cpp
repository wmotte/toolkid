#include "tkdHRFFit.h"
#include "tkdCmdParser.h"
#include <sstream>

int main( int argc, char ** argv )
{
  std::string inputFileName;
  std::string designFileName;
  std::string outputFileName;
  std::string outputFitFileName;
  std::string logFileName;
  std::vector< double > initial;
  std::vector< int > range;
  std::vector< int > baseline;
  bool fitAverage = false;
  bool relative = false;

  range.push_back( 0 );
  range.push_back( 25 );

  baseline.push_back( -10 );
  baseline.push_back( 0 );

  for( int i = 0; i < 3; ++i )
    {
    initial.push_back( 1 );
    }

  std::stringstream ss;
  ss << "Fit Gamma/HRF to fMRI data" << std::endl << std::endl;
  ss << "The fitted function depending on m, k, a, and b with parameter x is:" << std::endl;
  ss << "  g(x,m,k,a,b)= m + k * (((x/b)^(a-1) * exp(-x/b)) / (b*gamma(a)))" << std::endl;
  ss << "Here, m is the baseline value, k is a scalar multiplier, " << std::endl;
  ss << "and a and b determine the shape of the gamma function" << std::endl;

  tkd::CmdParser p( "hrffit",  ss.str() );

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "Input 1D vector text file" )
    ->SetRequired( true );

  p.AddArgument( designFileName, "design" )
    ->AddAlias( "d" )
    ->SetDescription( "Input design 1D vector text file (per volume stimulation off=0, on=1)" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "Output parameters (mean, start, end, k, alpha, beta)" )
    ->SetRequired( true );

  p.AddArgument( outputFitFileName, "output-fit" )
    ->AddAlias( "of" )
    ->SetDescription( "Output fitted vector" );

  p.AddArgument( logFileName, "log" )
    ->AddAlias( "l" )
    ->SetDescription( "Log file" );

  p.AddArgument( initial, "initial" )
    ->AddAlias( "n" )
    ->SetDescription( "Initial parameters: [ k alpha beta ]" )
    ->SetMinMax( 3, 3 );

  p.AddArgument( range, "range" )
    ->AddAlias( "r" )
    ->SetDescription( "Activation range [start, end). Default: [0, 25)" )
    ->SetMinMax( 2, 2 );

  p.AddArgument( baseline, "baseline" )
    ->AddAlias( "b" )
    ->SetDescription( "Baseline range [start, end). Default: [-10, 0)" )
    ->SetMinMax( 2, 2 );

  p.AddArgument( fitAverage, "average" )
    ->AddAlias( "a" )
    ->SetDescription( "Fit function to average stimulation" );

  p.AddArgument( relative, "relative" )
    ->AddAlias( "rel" )
    ->SetDescription( "Fit to data relative to baseline" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  tkd::HRFFit::VectorType parameters( initial.size() );
  for( int i = 0; i < initial.size(); ++i )
    {
    parameters( i ) = initial[ i ];
    }

  tkd::HRFFit fit;
  fit.SetInputFileName( inputFileName );
  fit.SetDesignFileName( designFileName );
  fit.SetOutputFileName( outputFileName );
  fit.SetOutputFitFileName( outputFitFileName );
  fit.SetLogFileName( logFileName );
  fit.SetBaseline( baseline[ 0 ], baseline[ 1 ] );
  fit.SetFitRange( range[ 0 ], range[ 1 ] );
  fit.SetInitialParameters( parameters );
  fit.SetFitAverage( fitAverage );
  fit.SetFitRelative( relative );
  fit.Run();

  return 0;
}
