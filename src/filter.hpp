#include "sparray.hpp"

#ifndef _DEMO_FILTER_H_
#define _DEMO_FILTER_H_

namespace demo {

template <class Pred>
sparray filter_seq(const Pred& p, sparray& xs) {
  long n = xs.size();
  long m = 0;
  for (long i = 0; i < n; i++) if (p(xs[i])) m++;

  sparray result(m);
  long j = 0;
  for (long i = 0; i < n; i++) {
    if (p(xs[i])) {
      result[j] = xs[i];
      j++;
    }
  }
  return result;
}

template <class Pred>
sparray filter(const Pred& p, const sparray& xs) {
  return sparray(0);
}

} /* namespace demo */

#endif /* _DEMO_FILTER_H_ */
