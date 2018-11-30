#include "benchmark.hpp"
#include "benchmark_type.hpp"
#include "sparray.hpp"
#include "fib.hpp"
#include "fib-ref.hpp"
#include "map.hpp"
#include "map-ref.hpp"
#include "filter.hpp"
#include "filter-ref.hpp"
#include "hash.hpp"
#include "sort.hpp"
#include <string>

benchmark_type fib_bench() {
  long n = pasl::util::cmdline::parse_or_default_long("n", 42);
  bool ref = pasl::util::cmdline::parse_or_default_bool("ref", false);
  long* result = new long;

  auto init = [=] { };

  auto bench = [=] {
#ifdef SEQUENTIAL_BASELINE
    *result = ref::fib_seq(n);
#else
    if (ref) *result = ref::fib_par(n);
    else     *result = demo::fib_par(n);
#endif
  };

  auto output = [=] {
    std::cout << "ref " << (ref ? "yes" : "no") << std::endl;
    std::cout << "result " << *result << std::endl;
  };

  auto destroy = [=] {
    delete result;
  };

  return make_benchmark(init, bench, output, destroy);
}

auto mapper_func = [] (long x) {
  return otherhash(otherhash(x));
};

benchmark_type map_bench() {
  long n = pasl::util::cmdline::parse_or_default_long("n", 1l<<20);
  bool ref = pasl::util::cmdline::parse_or_default_bool("ref", false);
  sparray* input = new sparray(0);
  long** result = new long*;

  auto init = [=] {
    *input = tabulate(identity_fct, n);
  };

  auto bench = [=] {
    long* data = &(*input)[0];
#ifdef SEQUENTIAL_BASELINE
    long* arr = new long[n];
    for (long i = 0; i < n; i++) arr[i] = mapper_func(data[i]);
    *result = arr;
#else
    if (ref) *result = ref::map(mapper_func, data, n);
    else     *result = demo::map(mapper_func, data, n);
#endif
  };

  auto output = [=] {
    std::cout << "result0 " << (*result)[0] << std::endl;
  };

  auto destroy = [=] {
    delete *result;
    delete result;
    delete input;
  };

  return make_benchmark(init, bench, output, destroy);
}

benchmark_type filter_bench() {
  long n = pasl::util::cmdline::parse_or_default_long("n", 1l<<20);
  bool ref = pasl::util::cmdline::parse_or_default_bool("ref", false);
  sparray* inp = new sparray(0);
  sparray* outp = new sparray(0);

  auto init = [=] {
    *inp = tabulate([] (long i) { return otherhash(i); }, n);
  };

  auto bench = [=] {
#ifdef SEQUENTIAL_BASELINE
    *outp = ref::filter_seq(is_even_fct, *inp);
#else
    if (ref) *outp = ref::filter(is_even_fct, *inp);
    else     *outp = demo::filter(is_even_fct, *inp);
#endif
  };

  auto output = [=] {
    std::cout << "result " << (*outp)[outp->size()-1] << std::endl;
  };

  auto destroy = [=] {
    //delete inp;
    //delete outp;
  };

  return make_benchmark(init, bench, output, destroy);
}

benchmark_type filter_sparray_bench() {
  long n = pasl::util::cmdline::parse_or_default_long("n", 1l<<20);
  sparray* inp = new sparray(0);
  sparray* outp = new sparray(0);

  auto init = [=] {
    *inp = tabulate([] (long i) { return otherhash(i); }, n);
  };

  auto bench = [=] {
#ifdef SEQUENTIAL_BASELINE
    *outp = ref::filter_seq(is_even_fct, *inp);
#else
    *outp = filter(is_even_fct, *inp);
#endif
  };

  auto output = [=] {
    std::cout << "result " << (*outp)[outp->size()-1] << std::endl;
  };

  auto destroy = [=] {
    //delete inp;
    //delete outp;
  };

  return make_benchmark(init, bench, output, destroy);
}

benchmark_type mergesort_bench() {
  long n = pasl::util::cmdline::parse_or_default_long("n", 1l<<20);
  sparray* inp = new sparray(0);
  sparray* outp = new sparray(0);

  auto init = [=] {
    *inp = tabulate([] (long i) { return otherhash(i); }, n);
  };

  auto bench = [=] {
#ifdef SEQUENTIAL_BASELINE
    *outp = seqsort(*inp);
#else
    *outp = mergesort_opt(*inp);
#endif
  };

  auto output = [=] {
    std::cout << "result " << (*outp)[0] << std::endl;
  };

  auto destroy = [=] {
    //delete inp;
    //delete outp;
  };

  return make_benchmark(init, bench, output, destroy);
}

/*---------------------------------------------------------------------*/
/* PASL Driver */

int main(int argc, char** argv) {

  benchmark_type bench;

  auto init = [&] {
    pasl::util::cmdline::argmap<std::function<benchmark_type()>> m;
    m.add("fib",                  [&] { return fib_bench(); });
    m.add("map",                  [&] { return map_bench(); });
    m.add("filter",               [&] { return filter_bench(); });
    m.add("filter-sparray",       [&] { return filter_sparray_bench(); });
    m.add("mergesort",            [&] { return mergesort_bench(); });
    bench = m.find_by_arg("bench")();
    bench_init(bench);
  };
  auto run = [&] (bool) {
    bench_run(bench);
  };
  auto output = [&] {
    bench_output(bench);
  };
  auto destroy = [&] {
    bench_destroy(bench);
  };
  pasl::sched::launch(argc, argv, init, run, output, destroy);
}

/***********************************************************************/
