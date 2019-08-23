#ifndef vnl_vector_to_std_vector_h_
#define vnl_vector_to_std_vector_h_

#include <vector>
#include <vnl/vnl_vector.h>

template <class T>
inline
std::vector<T> vnl_vector_to_std_vector(vnl_vector<T> const& cl)
{
  std::vector<T> l; l.reserve(cl.size());
  for (unsigned int i=0; i < cl.size(); ++i)
    l.push_back(cl[i]);
  return l;
}

#endif // vnl_vector_to_std_vector_h_
