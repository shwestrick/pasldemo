/* COPYRIGHT (c) 2018 Sam Westrick
 */

#include <math.h>
#include <climits>
#include "sparray.hpp"

#ifndef _MAP_DEMO_H_
#define _MAP_DEMO_H_

namespace demo {

template <class Func>
long* map(const Func& f, long* arr, long n) {
  long* result = new long[n];
  return result;
}

} /* namespace demo */

#endif /*! _MAP_DEMO_H_ */
