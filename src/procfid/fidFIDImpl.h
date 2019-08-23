#ifndef __fidFIDImpl_h__
#define __fidFIDImpl_h__
#define LINUX 1
#include "data.h"
#include "procparser.h"

namespace fid
{

class FIDImpl
{
public:
  datafilehead Header;
  Procparser Procpar;
};

} // end namespace fid

#endif /*__fidFIDImpl_h__*/
