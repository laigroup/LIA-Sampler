#ifndef _MERSENNE_CC
#define _MERSENNE_CC

#include "mersenne.h"

#include <ctime>

const int M = 397;
const unsigned int MATRIX_A = 0x9908b0dfUL;
const unsigned int UPPER_MASK = 0x80000000UL;
const unsigned int LOWER_MASK = 0x7fffffffUL;

Mersenne::Mersenne() {
    seed((int)std::time(0));
}

Mersenne::Mersenne(int s) {
    seed(s);
}

Mersenne::Mersenne(unsigned int *init_key, int key_length) {
    seed(init_key, key_length);
}

Mersenne::Mersenne(const Mersenne &copy) {
    *this = copy;
}

Mersenne &Mersenne::operator=(const Mersenne &copy) {
    for (int i = 0; i < N; i++)
        mt[i] = copy.mt[i];
    mti = copy.mti;
    return *this;
}

void Mersenne::seed(int se) {
    unsigned int s = ((unsigned int)(se << 1)) + 1;
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < N; mti++) {
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }
}

void Mersenne::seed(unsigned int *init_key, int key_length) {
    int i = 1, j = 0, k = (N > key_length ? N : key_length);
    seed(19650218UL);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) + init_key[j] + j;
        mt[i] &= 0xffffffffUL;
        i++;
        j++;
        if (i >= N) {
            mt[0] = mt[N - 1];
            i = 1;
        }
        if (j >= key_length)
            j = 0;
    }
    for (k = N - 1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i;
        mt[i] &= 0xffffffffUL;
        i++;
        if (i >= N) {
            mt[0] = mt[N - 1];
            i = 1;
        }
    }
    mt[0] = 0x80000000UL;
}

unsigned int Mersenne::next32() {
    unsigned int y;
    static unsigned int mag01[2] = {0x0UL, MATRIX_A};
    if (mti >= N) {
        int kk;
        for (kk = 0; kk < N - M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        mti = 0;
    }
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}

int Mersenne::next31() {
    return (int)(next32() >> 1);
}

double Mersenne::nextClosed() {
    unsigned int a = next32() >> 5, b = next32() >> 6;
    return (a * 67108864.0 + b) * (1.0 / 9007199254740991.0);
}

double Mersenne::nextHalfOpen() {
    unsigned int a = next32() >> 5, b = next32() >> 6;
    return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}

double Mersenne::nextOpen() {
    unsigned int a = next32() >> 5, b = next32() >> 6;
    return (0.5 + a * 67108864.0 + b) * (1.0 / 9007199254740991.0);
}

int Mersenne::next(int bound) {
    unsigned int value;
    do {
        value = next31();
    } while (value + (unsigned int)bound >= 0x80000000UL);
    return (int)(value % bound);
}

__int128_t Mersenne::random_num_in_range(__int128_t low, __int128_t high) {
    if (low > high) {
        std::swap(low, high);  // 确保 low <= high
    }

    unsigned int range = high - low + 1;
    unsigned int value;

    do {
        value = next31();
    } while (value >= 0x80000000UL - (0x80000000UL % range));

    return low + (value % range);
}

#endif
