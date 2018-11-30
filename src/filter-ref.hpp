/* COPYRIGHT (c) 2018 Sam Westrick
 */

#include <math.h>
#include <climits>
#include "sparray.hpp"

#ifndef _FILTER_REF_H_
#define _FILTER_REF_H_

namespace ref {

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
struct filter_contr {
  static loop_controller_type c;
};
template <class Pred>
loop_controller_type filter_contr<Pred>::c("filter_contr_" + par::string_of_template_arg<Pred>());

template <class Pred>
sparray filter(const Pred& p, const sparray& xs) {
  long n = xs.size();

  scan_excl_result offsets = scan_excl(
    [] (long a, long b) { return a + b; },
    p, // essentially a map-scan, with this function applied first
    0,
    xs);
  long m = offsets.total;

  sparray result(m);
  par::parallel_for(filter_contr<Pred>::c, (long)0, n, [&] (long i) {
    auto x = xs[i];
    long curr = offsets.partials[i];
    long next = i == n-1 ? m : offsets.partials[i+1];
    if (curr != next) {
      result[curr] = x;
    }
  });
  return result;
}

} /* namespace ref */

#endif /*! _FILTER_REF_H_ */
