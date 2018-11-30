/* COPYRIGHT (c) 2018 Sam Westrick
 */

#include <math.h>
#include <climits>
#include "sparray.hpp"

#ifndef _MAP_REF_H_
#define _MAP_REF_H_

namespace ref {

template <class Func>
struct map_contr {
  static loop_controller_type c;
};
template <class Func>
loop_controller_type map_contr<Func>::c("ref::map_contr_" + par::string_of_template_arg<Func>());

template <class Func>
long* map(const Func& f, long* arr, long n) {
  long* result = new long[n];
  par::parallel_for(map_contr<Func>::c, (long)0, n, [&] (long i) {
    result[i] = f(arr[i]);
  });
  return result;
}

} /* namespace ref */

#endif /*! _MAP_REF_H_ */
