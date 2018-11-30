/* COPYRIGHT (c) 2018 Sam Westrick
 */

#include <math.h>
#include <climits>
#include "sparray.hpp"

#ifndef _FIB_REF_H_
#define _FIB_REF_H_

namespace ref {

long fib_seq(long n) {
  if (n < 2) return n;
  return fib_seq(n-1) + fib_seq(n-2);
}

controller_type fib_par_contr("fib_par");

long fib_par(long n) {
  long result;
  par::cstmt(fib_par_contr, [&] { return 1l<<n; }, [&] {
    if (n < 2) result = n;
    else {
      long a,b;
      par::fork2([&] {
        a = fib_par(n-1);
      }, [&] {
        b = fib_par(n-2);
      });
      result = a+b;
    }
  }, [&] {
    result = fib_seq(n);
  });
  return result;
}

} /* namespace ref */

#endif /*! _FIB_REF_H_ */
