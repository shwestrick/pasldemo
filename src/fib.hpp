/* COPYRIGHT (c) 2018 Sam Westrick
 */

#include <math.h>
#include <climits>
#include "sparray.hpp"

#ifndef _FIB_DEMO_H_
#define _FIB_DEMO_H_

namespace demo {

long fib_seq(long n) {
  if (n < 2) return n;
  return fib_seq(n-1) + fib_seq(n-2);
}

long fib_par(long n) {
  return 0;
}

} /* namespace demo */

#endif /*! _FIB_DEMO_H_ */
