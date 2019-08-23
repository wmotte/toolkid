#ifndef __tkdHRFFit_h__
#define __tkdHRFFit_h__
#include "vnl/vnl_vector.h"

namespace tkd
{

class HRFFit
{
public:
  typedef double ScalarType;
  typedef vnl_vector< ScalarType > VectorType;

  HRFFit();

  /**
   * Perform fitting on average stimulation (default: false)
   */
  void SetFitAverage( bool fitAverage );

  /**
   * Perform fitting on data relative to the baseline (default: false)
   */
  void SetFitRelative( bool fitRelative );

  /**
   * 1D vector file, with 0 (stimulation off) or 1 (stimulation on) for each fMRI volume
   */
  void SetDesignFileName( const std::string& filename );

  /**
   * Log file
   */
  void SetLogFileName( const std::string& filename );

  /**
   * 1D vector file with fMRI data
   */
  void SetInputFileName( const std::string& filename );

  /**
   * Output parameters filename
   */
  void SetOutputFileName( const std::string& filename );

  /**
   * Output fitted vector filename
   */
  void SetOutputFitFileName( const std::string& filename );

  /**
   * Relative to first stimulation volume
   */
  void SetBaseline( int start, int end );

  /**
   * Relative to first stimulation volume
   */
  void SetFitRange( int start, int end );

  /**
   * 3 parameters: k, alpha, beta
   */
  void SetInitialParameters( const VectorType& parameters );

  /**
   * Perform fitting and save results
   */
  void Run();

protected:
  VectorType Read( const std::string& filename );

  std::string m_InputFileName;
  std::string m_OutputFileName;
  std::string m_OutputFitFileName;
  std::string m_DesignFileName;
  std::string m_LogFileName;

  int m_BaselineStart;
  int m_BaselineEnd;
  int m_FitStart;
  int m_FitEnd;

  bool m_FitAverage;
  bool m_FitRelative;

  VectorType m_Initial;
};

} // end namespace tkd

#endif /*__tkdHRFFit_h__*/
