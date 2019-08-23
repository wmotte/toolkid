#include "tkdHRFFit.h"
#include "tkdStimulus.h"
#include "tkdFitGammaHRF.h"
#include "tkdGammaHRF.h"
#include <iostream>
#include <fstream>
#include <vector>

namespace tkd
{

HRFFit::HRFFit() : m_FitAverage( false ), m_FitRelative( false )
{

}

void HRFFit::SetFitAverage( bool fitAverage )
{
  m_FitAverage = fitAverage;
}

void HRFFit::SetFitRelative( bool fitRelative )
{
  m_FitRelative = fitRelative;
}

void HRFFit::SetLogFileName( const std::string& filename )
{
  m_LogFileName = filename;
}

void HRFFit::SetDesignFileName( const std::string& filename )
{
  m_DesignFileName = filename;
}

void HRFFit::SetInputFileName( const std::string& filename )
{
  m_InputFileName = filename;
}

void HRFFit::SetOutputFileName( const std::string& filename )
{
  m_OutputFileName = filename;
}

void HRFFit::SetOutputFitFileName( const std::string& filename )
{
  m_OutputFitFileName = filename;
}

void HRFFit::SetBaseline( int start, int end )
{
  m_BaselineStart = start;
  m_BaselineEnd = end;
}

void HRFFit::SetFitRange( int start, int end )
{
  m_FitStart = start;
  m_FitEnd = end;
}

void HRFFit::SetInitialParameters( const VectorType& parameters )
{
  m_Initial = parameters;
}

HRFFit::VectorType HRFFit::Read( const std::string& filename )
{
  std::vector< ScalarType > data;
  std::ifstream in( filename.c_str() );
  while( !in.eof() && in.good() )
    {
    ScalarType x;
    in >> x;
    if ( in.good() )
      {
      data.push_back( x );
      }
    }

  VectorType vector = VectorType( data.size() );
  for( unsigned int i = 0; i < data.size(); ++i )
    {
    vector( i ) = data[ i ];
    }

  return vector;
}

void HRFFit::Run()
{
  VectorType data = Read( m_InputFileName );
  VectorType design = Read( m_DesignFileName );

  if ( design.size() > data.size() )
    {
    design = design.extract( data.size(), 0 );
    }

  std::ofstream logFile;
  std::ostream* log = &( std::cout );

  if ( m_LogFileName != "" )
    {
    logFile.open( m_LogFileName.c_str() );
    log = &( logFile );
    }

  std::vector< int > starts;
  std::vector< int > ends;
  std::vector< int > counts;

  for( int i = 0; i < design.size(); ++i )
    {
    if ( design( i ) > 0 && ( i == 0 || design( i - 1 ) == 0 ) )
      {
      starts.push_back( i );
      ends.push_back( i + 1 );
      counts.push_back( 1 );
      }
    else if ( ends.size() > 0 )
      {
      if ( design( i ) > 0 && design( i - 1 ) > 0 )
        {
        counts[ counts.size() - 1 ]++;
        }

      ends[ ends.size() - 1 ] = i + 1;
      }
    }

  std::vector< Stimulus > stimuli;
  for( int i = 0; i < starts.size(); ++i )
    {
    Stimulus stimulus;
    stimulus.BaselineStart = starts[ i ] + m_BaselineStart;
    stimulus.BaselineEnd = starts[ i ] + m_BaselineEnd;
    stimulus.StimulationStart = starts[ i ];
    stimulus.StimulationEnd = counts[ i ];
    stimulus.BlockStart = stimulus.BaselineStart;
    stimulus.BlockEnd = ends[ i ];
    stimulus.FitStart = starts[ i ] + m_FitStart;
    stimulus.FitEnd = starts[ i ] + m_FitEnd;

    if ( i == 0 )
      {
      stimulus.BlockStart = 0;
      }

    if ( ( i + 1 ) < starts.size() )
      {
      stimulus.BlockEnd += m_BaselineStart;
      }
    else
      {
      stimulus.BlockEnd = data.size();
      }

    stimuli.push_back( stimulus );
    }

  if ( m_FitAverage )
    {
    int minBefore, minAfter;

    for( int i = 0; i < stimuli.size(); ++i )
      {
      Stimulus stimulus = stimuli[ i ];
      int before = stimulus.BaselineEnd - stimulus.BlockStart;
      int after = stimulus.BlockEnd - stimulus.BaselineEnd;

      minBefore = ( ( i == 0 || minBefore > before ) ? before : minBefore );
      minAfter  = ( ( i == 0 || minAfter  > after  ) ? after  : minAfter );
      }

    VectorType average;
    ScalarType averageMean = 0;
    for( int i = 0; i < stimuli.size(); ++i )
      {
      Stimulus stimulus = stimuli[ i ];

      int from = stimulus.BaselineEnd - minBefore;
      int to = stimulus.BaselineEnd + minAfter;

      VectorType baseline = data.extract( stimulus.BaselineEnd - stimulus.BaselineStart, stimulus.BaselineStart );
      VectorType extract = data.extract( to - from, from );

      ScalarType mean = 0;
      for( int j = 0; j < baseline.size(); ++j )
        {
        mean += baseline( j );
        }
      mean /= static_cast< ScalarType >( baseline.size() );

      if ( average.size() == 0 )
        {
        average = ( extract - mean );
        }
      else
        {
        average += ( extract - mean );
        }

      averageMean += mean;
      }
    average /= static_cast< ScalarType >( stimuli.size() );
    averageMean /= static_cast< ScalarType >( stimuli.size() );
    average += averageMean;

    Stimulus stimulus;
    stimulus.BlockStart = 0;
    stimulus.BlockEnd = average.size();
    stimulus.BaselineStart = 0;
    stimulus.BaselineEnd = minBefore;
    stimulus.FitStart = minBefore + m_FitStart;
    stimulus.FitEnd = minBefore + m_FitEnd;
    stimulus.StimulationStart = stimulus.FitStart;
    stimulus.StimulationEnd = stimulus.FitStart;

    stimuli.clear();
    stimuli.push_back( stimulus );
    data = average;
    }

  std::ofstream out( m_OutputFileName.c_str() );
  out << "\"Stim\",\"k\",\"alpha\",\"beta\",\"mean\",\"SSD\",\"TTP\",\"Max\",\"X2\",\"P\",\"AUC_rel\",\"Max_rel\"" << std::endl;

  std::ofstream outFit;

  if ( m_OutputFitFileName != "" )
    {
    outFit.open( m_OutputFitFileName.c_str() );
    }

  VectorType fitted = VectorType( data.size() );
  fitted.fill( 0 );
  for( int i = 0; i < stimuli.size(); ++i )
    {
    Stimulus stimulus = stimuli[ i ];
    double mean = 0;
    int counter = 0;

    for( int j = stimulus.BaselineStart; j < stimulus.BaselineEnd; ++j, ++counter )
      {
      mean += data( j );
      }
    mean /= static_cast< double >( counter );

    ( *log ) << "# Block " << i << ": [" << stimulus.BlockStart << ", " << stimulus.BlockEnd << ") baseline [" << stimulus.BaselineStart << ", " << stimulus.BaselineEnd << "), fit [" << stimulus.FitStart << ", " << stimulus.FitEnd << "), mean: " << mean << std::endl;

    VectorType currentData = data.extract( stimulus.BlockEnd - stimulus.BlockStart, stimulus.BlockStart );
    VectorType currentParameters = m_Initial;
    VectorType currentFit = VectorType( currentData.size() );

    if ( m_FitRelative )
      {
      currentData = currentData / mean - 1.;
      mean = 0;
      }

    FitGammaHRF fit( mean,
                     stimulus.FitStart - stimulus.BlockStart,
                     stimulus.FitEnd - stimulus.BlockStart );

    ScalarType error = fit.Fit( currentData, currentParameters, currentFit, *log );

    GammaHRF function( currentData, mean, stimulus.FitStart - stimulus.BlockStart, currentData.size() );
    ScalarType timeToPeak = function.TimeToPeak( currentParameters );
    ScalarType maximumValue = function.Evaluate( timeToPeak, currentParameters );
    ScalarType chiSquared = function.ChiSquared( currentParameters );
    ScalarType chiSquaredP = function.ChiSquaredCDF( chiSquared );

    out << i << "," << currentParameters( 0 )
             << "," << currentParameters( 1 )
             << "," << currentParameters( 2 )
             << "," << mean
             << "," << error
             << "," << timeToPeak
             << "," << maximumValue
             << "," << chiSquared
             << "," << chiSquaredP
             << "," << ( currentParameters( 0 ) / mean )
             << "," << ( maximumValue / mean - 1. )
             << std::endl;

    for( int j = stimulus.BlockStart; j < stimulus.BlockEnd; ++j )
      {
      if ( m_FitRelative )
        {
        data( j ) = currentData( j - stimulus.BlockStart );
        }

      fitted( j ) = currentFit( j - stimulus.BlockStart );
      }
    }

  if ( m_OutputFitFileName != "" )
    {
    outFit << "\"frame\",\"data\",\"fit\"" << std::endl;
    for( int j = stimuli[ 0 ].BlockStart; j < stimuli[ stimuli.size() - 1 ].BlockEnd; ++j )
      {
      outFit << j << "," << data( j ) << "," << fitted( j ) << std::endl;
      }
    }
}

} // end namespace tkd
