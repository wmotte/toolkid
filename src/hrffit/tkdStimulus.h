#ifndef __tkdStimulus_h__
#define __tkdStimulus_h__

namespace tkd
{

class Stimulus
{
public:
  /**
   * Part of the data before stimulation, to calculate baseline value
   */
  int BaselineStart;
  int BaselineEnd;

  /**
   * Stimulation volumes
   */
  int StimulationStart;
  int StimulationEnd;

  /**
   * Part of the data to restrict the fitting to
   */
  int FitStart;
  int FitEnd;

  /**
   * Part of the data including baseline, stimulation, and response
   */
  int BlockStart;
  int BlockEnd;
};

} // end namespace tkd

#endif /*__tkdStimulus_h__*/
