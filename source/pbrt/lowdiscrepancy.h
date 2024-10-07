#pragma once
#include <vector>
#include <memory>
#include "pbrt/primes.h"
#include "pbrt/hash.h"

inline int PermutationElement(uint32_t i, uint32_t l, uint32_t p) {
    uint32_t w = l - 1;
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;
    do {
        i ^= p;
        i *= 0xe170893d;
        i ^= p >> 16;
        i ^= (i & w) >> 4;
        i ^= p >> 8;
        i *= 0x0929eb3f;
        i ^= p >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | p >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;
    } while (i >= l);
    return (i + p) % l;
}

class CDigitPermutation
{
public:
    CDigitPermutation() = default;

    CDigitPermutation(int base, uint32_t seed) :base(base)
    {
        n_digits = 0;
        float invBase = 1.0 / float(base);
        float inv_base_m = 1;
        while (1 - (base - 1) * inv_base_m < 1) {
            ++n_digits;
            inv_base_m *= invBase;
        }

        permutations.resize(n_digits * base);
        // Compute random permutations for all digits
        for (int digitIndex = 0; digitIndex < n_digits; ++digitIndex) {
            uint64_t dseed = pbrt::Hash(base, digitIndex, seed);
            for (int digitValue = 0; digitValue < base; ++digitValue) {
                int index = digitIndex * base + digitValue;
                permutations[index] = PermutationElement(digitValue, base, dseed);
            }
        }
    }

    int Permute(int digitIndex, int digitValue) const {
        return permutations[digitIndex * base + digitValue];
    }

private:
    int base, n_digits;
    std::vector<uint16_t> permutations;
};

using CDigitPermutationArrayPtr = std::shared_ptr<std::vector<CDigitPermutation>>;

CDigitPermutationArrayPtr computeRadicalInversePermutation(uint32_t seed)
{
	CDigitPermutationArrayPtr permutations = std::make_shared<std::vector<CDigitPermutation>>();
    permutations->resize(pbrt::PrimeTableSize);
    for (int i = 0; i < pbrt::PrimeTableSize; ++i)
        (*permutations)[i] = CDigitPermutation(pbrt::Primes[i], seed);
    return permutations;
}

inline float scrambledRadicalInverse(int base_idx,uint64_t a, const CDigitPermutation& permutation)
{
	uint32_t base = pbrt::Primes[base_idx];

	uint64_t limit = ~0ull / base - base;
	float inv_base = 1.0 / float(base);
	float inv_base_m = 1.0;
	uint64_t reverse_digits = 0;

	int digit_index = 0;
	while (1 - (base - 1) * inv_base_m < 1 && reverse_digits < limit)
	{
		uint64_t next = a / base;
		int digit_value = a - next * base;
		reverse_digits = reverse_digits * base + permutation.Permute(digit_index, digit_value);
		inv_base_m *= inv_base;
		++digit_index;
		a = next;
	}
	return std::min(inv_base_m * reverse_digits, float(0x1.fffffep-1));
}