// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_UTIL_PRIMES_H
#define PBRT_UTIL_PRIMES_H

// A7x[BEGIN]
//#include <pbrt/pbrt.h>
//#include <pbrt/util/pstd.h>
// A7x[END]

namespace pbrt {

// Prime Table Declarations
static constexpr int PrimeTableSize = 1000;
extern int Primes[PrimeTableSize];

}  // namespace pbrt

#endif  // PBRT_UTIL_PRIMES_H
