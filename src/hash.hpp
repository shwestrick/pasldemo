/* COPYRIGHT (c) 2014 Umut Acar, Arthur Chargueraud, and Michael Rainey
 * COPYRIGHT (c) 2018 Sam Westrick
 * All rights reserved.
 *
 * \file hash.hpp
 * \brief Hash function
 *
 */


#ifndef _MINICOURSE_HASH_H_
#define _MINICOURSE_HASH_H_

/***********************************************************************/

long otherhash(long i) {
  uint64_t v = ((uint64_t) i) * 3935559000370003845 + 2691343689449507681;
  v = v ^ (v >> 21);
  v = v ^ (v << 37);
  v = v ^ (v >> 4);
  v = v * 4768777513237032717;
  v = v ^ (v << 20);
  v = v ^ (v >> 41);
  v = v ^ (v <<  5);
  return (long) (v & ((((uint64_t) 1) << 63) - 1));
}

// taken from https://gist.github.com/badboy/6267743

long hash64shift(long key) {
  auto unsigned_right_shift = [] (long x, int y) {
    unsigned long r = (unsigned long) x >> y;
    return (long) r;
  };
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (unsigned_right_shift(key, 24));
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (unsigned_right_shift(key, 14));
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (unsigned_right_shift(key, 28));
  key = key + (key << 31);
  return key;
}

long hash_signed(long key) {
  return hash64shift(key);
}

unsigned long hash_unsigned(long key) {
  return (unsigned long)hash_signed(key);
}

unsigned long random_index(long key, long n) {
  unsigned long x = (unsigned long)hash64shift(key);
  return x % n;
}

int log2_up(unsigned long i) {
  int a=0;
  long b=i-1;
  while (b > 0) {b = b >> 1; a++;}
  return a;
}

/***********************************************************************/

#endif /*! _MINICOURSE_HASH_H_ */
