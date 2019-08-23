#ifndef std_vector_to_vnl_vector_h_
#define std_vector_to_vnl_vector_h_

#include <vector>
#include <vnl/vnl_vector.h>

template <class T>
inline
vnl_vector<T> std_vector_to_vnl_vector( std::vector<T> const& cl)
{
  vnl_vector<T> l( cl.size() );
  for (unsigned int i=0; i < cl.size(); ++i)
    l( i ) = cl.at( i );
  return l;
}

#endif // std_vector_to_vnl_vector_h_
