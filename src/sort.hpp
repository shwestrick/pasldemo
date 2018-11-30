/* COPYRIGHT (c) 2014 Umut Acar, Arthur Chargueraud, and Michael
 * Rainey
 * COPYRIGHT (c) 2018 Sam Westrick
 * All rights reserved.
 *
 * \file sort.hpp
 * \brief Sorting algorithms
 *
 */

#include <cstring>
#include <cmath>

//include "native.hpp"
#include "sparray.hpp"
//include "sparray_triples.hpp"

#ifndef _MINICOURSE_SORT_H_
#define _MINICOURSE_SORT_H_

/***********************************************************************/

long nlogn(long n) {
  return pasl::data::estimator::annotation::nlgn(n);
}

/*---------------------------------------------------------------------*/
/* Sequential sort */

void in_place_sort(sparray& xs, long lo, long hi) {
  long n = hi-lo;
  if (n < 2)
    return;
  std::sort(&xs[lo], &xs[hi-1]+1);
}

void in_place_sort(sparray& xs) {
  in_place_sort(xs, 0l, xs.size());
}

sparray seqsort(const sparray& xs) {
  sparray tmp = copy(xs);
  in_place_sort(tmp);
  return tmp;
}

sparray seqsort(const sparray& xs, long lo, long hi) {
  sparray tmp = slice(xs, lo, hi);
  in_place_sort(tmp);
  return tmp;
}

/*---------------------------------------------------------------------*/
/* Parallel quicksort */

value_type median(value_type a, value_type b, value_type c) {
  return  a < b ? (b < c ? b : (a < c ? c : a))
  : (a < c ? a : (b < c ? c : b));
}

controller_type quicksort_contr("quicksort");

sparray quicksort(const sparray& xs) {
  long n = xs.size();
  sparray result = { };
  auto seq = [&] {
    result = seqsort(xs);
  };
  par::cstmt(quicksort_contr, [&] { return nlogn(n); }, [&] {
    if (n <= 4) {
      seq();
    } else {
      value_type p = median(xs[n/4], xs[n/2], xs[3*n/4]);
      sparray less = filter([&] (value_type x) { return x < p; }, xs);
      sparray equal = filter([&] (value_type x) { return x == p; }, xs);
      sparray greater = filter([&] (value_type x) { return x > p; }, xs);
      sparray left = { };
      sparray right = { };
      par::fork2([&] {
        left = quicksort(less);
      }, [&] {
        right = quicksort(greater);
      });
      result = concat(left, equal, right);
    }
  }, seq);
  return result;
}

controller_type in_place_quicksort_contr("in_place_quicksort");

void in_place_quicksort_rec(value_type* A, long n) {
  if (n < 2) {
    return;
  }
  par::cstmt(in_place_quicksort_contr, [&] { return nlogn(n); }, [&] {
    value_type p = A[0];
    value_type* L = A;   // below L are less than pivot
    value_type* M = A;   // between L and M are equal to pivot
    value_type* R = A+n-1; // above R are greater than pivot
    while (true) {
      while (! (p < *M)) {
        if (*M < p) std::swap(*M,*(L++));
        if (M >= R) break;
        M++;
      }
      while (p < *R) R--;
      if (M >= R) break;
      std::swap(*M,*R--);
      if (*M < p) std::swap(*M,*(L++));
      M++;
    }
    par::fork2([&] {
      in_place_quicksort_rec(A, L-A);
    }, [&] {
      in_place_quicksort_rec(M, A+n-M); // Exclude all elts that equal pivot
    });
  }, [&] {
    std::sort(A, A+n);
  });
}

sparray in_place_quicksort(const sparray& xs) {
  sparray result = copy(xs);
  long n = xs.size();
  if (n == 0) {
    return result;
  }
  in_place_quicksort_rec(&result[0], n);
  return result;
}

// sparray quicksort_seqpart(const sparray& xs) {
//   long n = xs.size();
//   sparray result = { };
//   auto seq = [&] {
//     result = seqsort(xs);
//   };
//   par::cstmt(quicksort_contr, [&] { return nlogn(n); }, [&] {
//     if (n <= 4) {
//       seq();
//     } else {
//       value_type p = median(xs[n/4], xs[n/2], xs[3*n/4]);
//       // sparray less = filter([&] (value_type x) { return x < p; }, xs);
//       // sparray equal = filter([&] (value_type x) { return x == p; }, xs);
//       // sparray greater = filter([&] (value_type x) { return x > p; }, xs);
//       sparray left = { };
//       sparray right = { };
//       par::fork2([&] {
//         left = quicksort(less);
//       }, [&] {
//         right = quicksort(greater);
//       });
//       result = concat(left, equal, right);
//     }
//   }, seq);
//   return result;
// }


// triple_type add_triples(triple_type t1, triple_type t2) {
//   return std::make_tuple<long,long,long>(std::get<0>(t1) + std::get<0>(t2),
//                                          std::get<1>(t1) + std::get<1>(t2),
//                                          std::get<2>(t1) + std::get<2>(t2));
// }
//
// controller_type quicksort_fast_contr("quicksort_fast");
// loop_controller_type quicksort_fast_filter_contr("quicksort_fast_filter");
//
// sparray quicksort_fast(const sparray& xs) {
//   long n = xs.size();
//   sparray result = { };
//   auto seq = [&] {
//     result = seqsort(xs);
//   };
//   par::cstmt(quicksort_fast_contr, [&] { return nlogn(n); }, [&] {
//     if (n <= 4) {
//       seq();
//     } else {
//       value_type p = median(xs[n/4], xs[n/2], xs[3*n/4]);
//
//       auto against_pivot = [&] (value_type x) {
//         if (x < p) return std::make_tuple<long,long,long>(1,0,0);
//         else if (x == p) return std::make_tuple<long,long,long>(0,1,0);
//         else return std::make_tuple<long,long,long>(0,0,1);
//       };
//       sparray_triples trips = tabulate_triples([&] (long i) { return against_pivot(xs[i]); }, xs.size());
//       triple_type zero_triple = std::make_tuple<long,long,long>(0,0,0);
//       scan_triples_excl_result offsets = scan_triples_excl(add_triples, zero_triple, trips);
//
//       sparray less = sparray(std::get<0>(offsets.total));
//       sparray equal = sparray(std::get<1>(offsets.total));
//       sparray greater = sparray(std::get<2>(offsets.total));
//       par::parallel_for(quicksort_fast_filter_contr, 0l, xs.size(), [&] (long i) {
//         auto tripi = offsets.partials[i];
//         long curr0 = std::get<0>(tripi);
//         long curr1 = std::get<1>(tripi);
//         long curr2 = std::get<2>(tripi);
//         long next0, next1;
//         if (i + 1 == n) {
//           next0 = std::get<0>(offsets.total);
//           next1 = std::get<1>(offsets.total);
//         } else {
//           auto tripnext = offsets.partials[i+1];
//           next0 = std::get<0>(tripnext);
//           next1 = std::get<1>(tripnext);
//         }
//
//         if (curr0+1 == next0)
//           less[curr0] = xs[i];
//         else if (curr1+1 == next1)
//           equal[curr1] = xs[i];
//         else
//           greater[curr2] = xs[i];
//       });
//
//       // sparray less = filter([&] (value_type x) { return x < p; }, xs);
//       // sparray equal = filter([&] (value_type x) { return x == p; }, xs);
//       // sparray greater = filter([&] (value_type x) { return x > p; }, xs);
//       sparray left = { };
//       sparray right = { };
//       par::fork2([&] {
//         left = quicksort(less);
//       }, [&] {
//         right = quicksort(greater);
//       });
//       result = concat(left, equal, right);
//     }
//   }, seq);
//   return result;
// }

/*---------------------------------------------------------------------*/
/* Parallel mergesort */

// template <bool use_parallel_merge=true>
// sparray mergesort(const sparray& xs) {
//   // placeholder for reference solution
//   return copy(xs);
// }

template <bool use_parallel_merge=true>
sparray cilksort(const sparray& xs) {
  // placeholder for reference solution
  return copy(xs);
}

/***********************************************************************/

void merge_sequential(const sparray& xs, sparray& result, long result_lo,
                      long Alo, long Ahi, long Blo, long Bhi) {
  long i = Alo;
  long j = Blo;
  long k = result_lo;
  while (i < Ahi && j < Bhi) {
    if (xs[i] < xs[j]) {
      result[k] = xs[i];
      i++;
    } else {
      result[k] = xs[j];
      j++;
    }
    k++;
  }
  for (; i < Ahi; i++) { result[k] = xs[i]; k++; }
  for (; j < Bhi; j++) { result[k] = xs[j]; k++; }
}

controller_type mergesort_contr("mergesort");
controller_type mergesort_seqmerge_contr("mergesort_seqmerge");
void mergesort_rec(sparray& result, sparray& temp, long lo, long hi);
void mergesort_seqmerge_rec(sparray& result, sparray& temp, long lo, long hi);
void merge(const sparray& xs, sparray& result, long result_lo,
           long Alo, long Ahi, long Blo, long Bhi);


/*@brief: Wrapper for the recursive mergesort*/
sparray mergesort(const sparray& input) {
  long len = input.size();
  // Result start off as input and gets more sorted
  sparray result = copy(input);
  // Allocated a temp upfront, so we don't allocate anything later.
  sparray temp = copy(input);
  mergesort_rec(result, temp, 0L, len);
  return result;
}

sparray mergesort_seqmerge(const sparray& input) {
  long len = input.size();
  // Result start off as input and gets more sorted
  sparray result = copy(input);
  // Allocated a temp upfront, so we don't allocate anything later.
  sparray temp = copy(input);
  mergesort_seqmerge_rec(result, temp, 0L, len);
  return result;
}

/*@brief: Mergesorts result in place, and uses temp
 *        to hold temporary changes.
 *@note: in_place_sort() is defined in sort.hpp
 *@param: lo is inclusive!
 *@param: hi is exclusive!
 */
void mergesort_rec(sparray& result, sparray& temp,
                   long lo, long hi) {
    long len = hi - lo;
    if(len < 2) return;
    par::cstmt(mergesort_contr, [&] { return nlogn(len); }, [&] {
      long mid = lo + len / 2;
      par::fork2([&] {
        mergesort_rec(result, temp, lo, mid);
      }, [&] {
        mergesort_rec(result, temp, mid, hi);
      });
      merge(result, temp, lo, lo, mid, mid, hi);
      // copy from temp back into result
      prim::pcopy(&temp[0], &result[0], lo, hi, lo);
    }, [&] {
      in_place_sort(result, lo, hi); // Leave this alone!
    });
}

void mergesort_seqmerge_rec(sparray& result, sparray& temp,
                   long lo, long hi) {
    long len = hi - lo;
    if(len < 2) return;
    par::cstmt(mergesort_seqmerge_contr, [&] { return nlogn(len); }, [&] {
      long mid = lo + len / 2;
      par::fork2([&] {
        mergesort_rec(result, temp, lo, mid);
      }, [&] {
        mergesort_rec(result, temp, mid, hi);
      });
      merge_sequential(result, temp, lo, lo, mid, mid, hi);
      // copy from temp back into result
      prim::pcopy(&temp[0], &result[0], lo, hi, lo);
    }, [&] {
      in_place_sort(result, lo, hi); // Leave this alone!
    });
}

void merge(const sparray& xs, sparray& result, long result_lo,
           long Alo, long Ahi, long Blo, long Bhi) {
  return;
}

/***********************************************************************/

void serial_put_sort(const sparray& src, sparray& dst, long lo, long hi) {
  prim::copy(&src[0], &dst[0], lo, hi, lo);
  in_place_sort(dst, lo, hi);
}

controller_type mergesort_dc_rec_contr("mergesort_dc_rec_contr");

// If stay, sorts A[lo,hi) in place. Otherwise, puts the result in B[lo,hi).
void mergesort_dc_rec(sparray& A, sparray& B, long lo, long hi, bool stay) {
  long n = hi - lo;
  auto seq = [&] {
    if (stay) in_place_sort(A, lo, hi);
    else serial_put_sort(A, B, lo, hi);
  };
  if (n < 2) seq();
  else par::cstmt(mergesort_dc_rec_contr, [&] { return nlogn(n); }, [&] {
    long mid = lo + n / 2;
    par::fork2([&] {
      mergesort_dc_rec(A, B, lo, mid, !stay);
    }, [&] {
      mergesort_dc_rec(A, B, mid, hi, !stay);
    });
    if (stay) merge(B, A, lo, lo, mid, mid, hi);
    else merge(A, B, lo, lo, mid, mid, hi);
  }, seq);
}

sparray mergesort_dc(const sparray& xs) {
  long n = xs.size();
  sparray A = copy(xs);
  sparray B = sparray(n);
  mergesort_dc_rec(A, B, 0, n, true);
  return A;
}

sparray mergesort_opt(const sparray& xs) {
  return mergesort_dc(xs);
}

/***********************************************************************/


// loop_controller_type copy_slice_contr("copy_slice");
// void copy_slice(sparray& from, sparray& to, long flo, long fhi, long tlo){
//     par::parallel_for(copy_slice_contr, 0L, fhi-flo, [&] (long i){
//             to[tlo + i] = from[flo + i];
//         });
// }

// void ref_serial_merge(sparray& input, sparray& temp, long start,
//                       long Alo, long Ahi, long Blo, long Bhi) {
void ref_serial_merge(const sparray& xs, sparray& result, long result_lo,
                      long Alo, long Ahi, long Blo, long Bhi) {
  long i = Alo;
  long j = Blo;
  long k = result_lo;
  while (i < Ahi && j < Bhi) {
    if (xs[i] < xs[j]) {
      result[k] = xs[i];
      i++;
    } else {
      result[k] = xs[j];
      j++;
    }
    k++;
  }
  for (; i < Ahi; i++) { result[k] = xs[i]; k++; }
  for (; j < Bhi; j++) { result[k] = xs[j]; k++; }
}

std::pair<long,long> ref_dual_kth_smallest(const sparray& input, long Alo, long Ahi,
                                           long Blo, long Bhi, long k){
    long Alen = Ahi - Alo;
    long Blen = Bhi - Blo;
    std::pair<long, long> res;
    assert(Alen + Blen >= k);
    assert(k >= 0);

    // Use A to specify exploration range
    // Thus use the smaller one as A
    if(Alen > Blen){
        res = ref_dual_kth_smallest(input, Blo, Bhi, Alo, Ahi, k);
        return std::pair<long, long>(res.second, res.first);
    }

    // This will cause A_idx = -1, B_idx = -1. Not okay
    if(k == 0){
        return std::pair<long, long>(0,0);
    }

    // Algorithm works by narrowing the range of possible A_idx values
    long A_idx, B_idx;
    long lower = ((long)k) - 1 - ((long)Blen);
    long lower_bound = std::max(-1L, lower); // Lower is inclusive
    long upper_bound = std::min(Alen, (long)k); // Upper is exclusive
    while(lower_bound < upper_bound){
        A_idx = ((upper_bound - lower_bound) / 2) + lower_bound;
        B_idx = k-2 - A_idx;
        // Case 1: Prove A_idx is too high
        if((A_idx != -1) && // -1 can't be too high
           ((B_idx == -1) || // A contains everyone
            (input[Alo + A_idx] > input[Blo + B_idx]))){

            if(B_idx != Blen - 1 && // We can't raise B_idx anymore
               input[Alo + A_idx] > input[Blo + B_idx + 1]){
                upper_bound = A_idx;
                continue;
            } else {
                return std::pair<long, long>(A_idx + 1, B_idx + 1);
            }

        } else { // Case 2: Prove A_idx is too low

            if(A_idx != Alen - 1 && // We can lower A_idx
               input[Blo + B_idx] > input[Alo + A_idx + 1]){
                lower_bound = A_idx + 1;
                continue;
            } else {
                return std::pair<long, long>(A_idx + 1, B_idx + 1);
            }

        }
    }
    return std::pair<long, long>(lower_bound + 1, k-1 - lower_bound);
}

loop_controller_type merge_logn_idx_contr("merge_logn_idx");
loop_controller_type merge_logn_contr("merge_logn");

void ref_merge_log(const sparray& input, sparray& result, long result_lo,
                   long Astart, long Aend, long Bstart, long Bend) {
  long Alen = Aend - Astart;
  long Blen = Bend - Bstart;

  // Break the array into chunks of size n/(log(n))
  long part_size = floor(log2((float)(Alen + Blen)));
  long nparts = ceil((float)(Alen + Blen)/(float)part_size);

  // For the ith block, the slice we want from A is A[i] to A[i+1]
  std::pair<long,long>* idxs = my_malloc<std::pair<long,long>>(nparts + 1);

  long search_size = floor(log2((float)(std::min(Alen, Blen))));
  auto compl_fct_idx = [&] (long lo, long hi) {
    return (hi - lo) * search_size;
  };
  par::parallel_for(merge_logn_idx_contr, compl_fct_idx, (long)0, nparts, [&] (long i) {
    idxs[i] = ref_dual_kth_smallest(input, Astart, Aend,
                                    Bstart, Bend, i * part_size);
  });
  idxs[nparts] = std::make_pair(Alen,Blen);

  //copy_slice(input, temp, 0L, input.size(), 0L);
  // Merge the two slices, and copy into the correct location
  auto compl_fct_merge = [&] (long lo, long hi) {
    return (hi - lo) * part_size;
  };
  par::parallel_for(merge_logn_contr, compl_fct_merge, (long)0, nparts, [&] (long j) {
    auto j_idxs = idxs[j];
    auto jnext_idxs = idxs[j+1];
//    long start = j_idxs.first + j_idxs.second + std::min(Astart, Bstart);
    ref_serial_merge(input, result, result_lo + j_idxs.first + j_idxs.second,
                     Astart + j_idxs.first, Astart + jnext_idxs.first,
                     Bstart + j_idxs.second, Bstart + jnext_idxs.second);
  });
  free(idxs);
  return;
}


controller_type mergesort_logn_contr("mergesort_logn");

// void ref_mergesort_inplace(sparray& input, sparray& temp, long start, long end){
//     long len = end - start;
//     if(len < 2) return;
//     par::cstmt(mergesort_logn_contr, [&] { return nlogn(len); }, [&] {
//       par::fork2([&] {
//               ref_mergesort_inplace(input, temp, start, start + len/2);
//           }, [&] {
//               ref_mergesort_inplace(input, temp, start + len/2, end);
//           });
//       ref_merge_log(input, temp, start, start + len/2, start + len/2, end);
//       prim::pcopy(&temp[0], &result[0], lo, hi, lo);
//     }, [&] {
//       in_place_sort(input, start, end);
//     });
//     return;
// }

void ref_mergesort_rec(sparray& result, sparray& temp,
                   long lo, long hi) {
  long len = hi - lo;
  if(len < 2) return;
  par::cstmt(mergesort_logn_contr, [&] { return nlogn(len); }, [&] {
    long mid = lo + len / 2;
    par::fork2([&] {
      ref_mergesort_rec(result, temp, lo, mid);
    }, [&] {
      ref_mergesort_rec(result, temp, mid, hi);
    });
    ref_merge_log(result, temp, lo, lo, mid, mid, hi);
    // copy from temp back into result
    prim::pcopy(&temp[0], &result[0], lo, hi, lo);
  }, [&] {
    in_place_sort(result, lo, hi); // Leave this alone!
  });
}

sparray mergesort_logn(const sparray& input){
  sparray res = copy(input);
  sparray temp = copy(input);
  ref_mergesort_rec(res, temp, (long)0, input.size());
  return res;
}

#endif /*! _MINICOURSE_SORT_H_ */
