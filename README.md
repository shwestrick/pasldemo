# PASL Demo

Below is a sketch for a walkthrough of using PASL
(https://github.com/deepsea-inria/pasl). These are basically my lecture
notes for giving a demonstration of PASL in the course 15-210 at CMU.

# Intro

```sh
make demo.opt
./demo.opt -bench fib --ref -n 42 -proc 1
./demo.opt -bench fib --ref -n 42 -proc 72

make demo.log
./demo.log -bench fib --ref -n 42 -proc 72 --pview
```

# Fibonacci

A simple parallel implementation:

```c++
long fib_par(long n) {
  if (n < 2) return n;
  long a,b;
  par::fork2([&] {
    a = fib_par(n-1);
  }, [&] {
    b = fib_par(n-2);
  });
  return a + b;
}
```

Let's run it.

```sh
make demo.opt
./demo.opt -bench fib -n 33 -proc 1
./demo.opt -bench fib -n 33 -proc 72

prun speedup -baseline "./demo.opt" -parallel "./demo.opt -proc 1,18,36,54,72" \
-bench fib -n 33
```

Okay, it's parallel, some speedup... we're done?

*How do we know it's performing well?* We need to compare against a
**sequential baseline**.

```c++
long fib_seq(long n) {
  if (n < 2) return n;
  return fib_seq(n-1) + fib_seq(n-2);
}
```

```sh
make demo.baseline
./demo.baseline -bench fib -n 33
```

Factor 30 observed work efficiency. It's slow!

We need **granularity control**. Here's an example using
a constant threshold ("manual" granularity control).

```c++
const long threshold = 20;
long fib_par(long n) {
  if (n <= threshold) return fib_seq(n);
  long a,b;
  par::fork2([&] {
    a = fib_par(n-1);
  }, [&] {
    b = fib_par(n-2);
  });
  return a + b;
}
```

```sh
make demo.opt
./demo.opt -bench fib -n 33
./demo.baseline -bench fib -n 33
./demo.opt -bench fib -n 42
./demo.baseline -bench fib -n 42

prun speedup -baseline "./demo.opt" -parallel "./demo.opt -proc 1,18,36,54,72" \
-bench fib -n 42
```

Much better. No significant overhead, and it scales better too!

Manual granularity control in small examples might be acceptable, but it
quickly becomes unwieldy and error-prone. In particular, constant thresholding
is not:
  * *portable* (must be tuned for each machine)
  * *composable* (particularly in examples of higher-order primitives)

PASL provides primitives for *automatic granularity control* with complexity
annotations. Here we use `par::cstmt` which takes a complexity function, a
parallel body, and an alternative sequential body. It also takes a "controller"
which is a piece of PASL magic. Each controller must be unique to a particular
callsite (we can't reuse controllers for other functions).

```c++
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
```

```sh
make demo.opt
prun speedup -baseline "./demo.opt" -parallel "./demo.opt -proc 1,18,36,54,72" \
-bench fib -n 42
```

# Map

Using what we've seen so far, let's write map. We
use C++ lambdas and templating for the higher-order behavior.

```c++
template <class Func>
void* map_range(const Func& f, long* arr, long* result, long lo, long hi) {
  auto seq = [&] {
    for (long i = lo; i < hi; i++) result[i] = f(arr[i]);
  };
  par::cstmt(map_contr<Func>::c, [&] { return hi - lo; }, [&] {
    if (hi - lo < 2) seq();
    else {
      long mid = lo + (hi - lo)/2;
      par::fork2([&] { map_range(f, arr, result, lo, mid); },
                 [&] { map_range(f, arr, result, mid, hi); });
    }
  }, seq);
}

template <class Func>
long* map(const Func& f, long* arr, long n) {
  long* result = new long[n];
  map_range(f, arr, result, 0, n);
  return result;
}
```

The controller can be digested simply as a magic incantation.

```c++
template <class Func>
struct map_contr {
  static controller_type c;
};
template <class Func>
controller_type map_contr<Func>::c("map_contr_" + par::string_of_template_arg<Func>());
```

It gets good performance (this example mapping a hash function at each element):

```sh
make demo.baseline demo.opt
prun speedup -baseline "./demo.baseline" -parallel \
"./demo.opt -proc 1,18,36,54,72" -bench map -n 200000000 -runs 2
```

## Parallel-For

The code for `map` exhibits a common pattern which can be summarized as an
instance of a "parallel for loop". PASL provides a primitive for this. The
following code is approximately equivalent to what we wrote before:

```c++
template <class Func>
long* map(const Func& f, long* arr, long n) {
  long* result = new long[n];
  par::parallel_for(map_contr<Func>::c, (long)0, n, [&] (long i) {
    result[i] = f(arr[i]);
  });
  return result;
}
```

And the controller is, once again, a magic incantation (slightly different
than before, with `loop_controller_type` instead of `controller_type`).

```c++
template <class Func>
struct map_contr {
  static loop_controller_type c;
};
template <class Func>
loop_controller_type map_contr<Func>::c("map_contr_" + par::string_of_template_arg<Func>());
```

We can repeat the experiment to see basically the same performance.

```sh
make demo.baseline demo.opt
prun speedup -baseline "./demo.baseline" -parallel \
"./demo.opt -proc 1,18,36,54,72" -bench map -n 200000000 -runs 2
```

# Simple Parallel Arrays

Using these ideas, we provide a simple library reminiscent of ArraySequences in
210. The arrays are not polymorphic. Here are some examples:

```c++
sparray a = sparray(n);
sparray b = tabulate([&] (long i) { return a[i] + 1; }, n);
sparray c = map([&] (value_type x) { return x + 1; }, b);
sparray d = concat(b, c);
long m = d.size();
value_type sum = reduce([] (value_type x, value_type y) { return x + y; }, 0L, s);
scan_excl_result r = scan_excl([] (value_type x, value_type y) { return x + y; }, 0L, s);
// r.partials : sparray
// r.total : value_type
```

# Exercise: Filter

First, sequentially.

```c++
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
```

How to do in parallel?

```c++
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
```

```c++
template <class Pred>
struct filter_contr {
  static loop_controller_type c;
};
template <class Pred>
loop_controller_type filter_contr<Pred>::c("filter_contr_" + par::string_of_template_arg<Pred>());
```

Experiments:

```sh
make demo.baseline demo.opt

prun speedup -baseline "./demo.baseline" -parallel \
"./demo.opt -proc 1,18,36,54,72" -bench filter -n 200000000 -runs 2

prun speedup -baseline "./demo.baseline" -parallel \
"./demo.opt -proc 1,18,36,54,72" -bench filter-sparray -n 200000000 -runs 2

prun speedup -baseline "./demo.baseline" -parallel \
"./demo.opt -proc 1,18,36,54,72" -bench filter,filter-sparray -n 200000000 -runs 2
pplot speedup -series bench

make demo.log
./demo.log -bench filter -n 200000000 -proc 72 --pview
pview
./demo.log -bench filter-sparray -n 200000000 -proc 72 --pview
pview
```
