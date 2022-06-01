#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "wave.h"

#define PRNG_SHA3

#ifdef PRNG_SHA3
void FIPS202_SHA3_512(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);
#endif

void fill(prng_t PRNG) {
	FIPS202_SHA3_512(PRNG->buff, 64, PRNG->buff);
	PRNG->available = 64;
}

void fill2(prng_t PRNG) {
	PRNG->buff2 = rnd8(PRNG);
	PRNG->available2 = 8;
}

void fill3(prng_t PRNG) {
	// If the random byte is in [243,255], the first 5 trits of its
	// decomposition in base 3 are not uniformly distributed. We discard
	// those bytes (5%)
	do {
		PRNG->buff3 = rnd8(PRNG);
	} while (PRNG->buff3 >= 243); // 243 = 3^5
	PRNG->available3 = 5;
}

uint8_t rnd8(prng_t PRNG) {
	PRNG->bytecount++;
#ifdef PRNG_SHA3
	if (PRNG->available == 0) {
		fill(PRNG);
	}
	PRNG->available--;
	return PRNG->buff[PRNG->available];
#else
	return random() & 0xff;
#endif
}

uint8_t rnd_bit(prng_t PRNG) {
	uint8_t c;
	if (PRNG->available2 == 0) {
		fill2(PRNG);
	}
	PRNG->available2--;
	c = PRNG->buff2 & 1;
	PRNG->buff2 >>= 1;
	return c;
}

uint8_t rnd_trit(prng_t PRNG) {
	uint8_t c;
	if (PRNG->available3 == 0) {
		fill3(PRNG);
	}
	PRNG->available3--;
	c = PRNG->buff3 % 3;
	PRNG->buff3 /= 3;
	return c;
}

uint16_t rnd16(prng_t PRNG) {
	uint16_t r = rnd8(PRNG);
	r <<= 8;
	r ^= rnd8(PRNG);
	return r;
}

/* Random integer uniform between 0 and n-1 included.
	 We assume n < 2^16 (else equivalent to rnd16())
 */
uint16_t rnd_short(int n, prng_t PRNG) {
	uint16_t r, max;
	/* (max + 1) = l * n, the largest multiple of n <= 2^16 */
	max = 0xffff - (0x10000 % n);
	/* accept any r in [0, max] */
	/* if n >= 2^16 then max = -1 = 2^16 - 1 */
	do {
		r = rnd16(PRNG);
	} while (r > max); // else r is not uniformly distributed
	return r % n;
}

prng_t prng_init(char * type, char * params_id, unsigned long seed) {
	prng_t PRNG = malloc(sizeof (struct prng));
#ifdef PRNG_SHA3
	// 20 is the decimal size of the largest 'unsigned long'
	PRNG->init_str = malloc(strlen(type) + strlen(params_id) + 20 + 1);
	sprintf(PRNG->init_str, "%s%s%lu", type, params_id, seed);
	FIPS202_SHA3_512((unsigned char *) PRNG->init_str, strlen(PRNG->init_str), PRNG->buff);
	PRNG->available = 64;
#else
	PRNG->init_str = NULL;
	srandom(seed);
#endif
	PRNG->available2 = PRNG->available3 = 0;
	PRNG->bytecount = 0;
	return PRNG;
}

void prng_clear(prng_t PRNG) {
	free(PRNG->init_str);
	free(PRNG);
}

/* randomly shuffle the array a[] of length n so that the first t
	 coordinates is a random uniform subset of size t of the array
	 values */
void rand_shuffle(int * a, int n, int t, prng_t PRNG) {
	int i, r, c;

	for (i = 0; i < t; ++i) {
		r = i + rnd_short(n - i, PRNG); /* i <= r < n uniform */
		c = a[r]; a[r] = a[i]; a[i] = c;
	}
}

int * identity_perm(int n) {
	int * a = (int *) malloc(n * sizeof (int));
	for (int i = 0; i < n; ++i) a[i] = i;
	return a;
}

int * randperm(int n, prng_t PRNG) {
	int * a = identity_perm(n);
	rand_shuffle(a, n, n, PRNG);
	return a;
}

int * randcw(int n, int t, prng_t PRNG) {
	int * a = identity_perm(n);
	rand_shuffle(a, n, t, PRNG);
	return a;
}

int * rand_bincw(int n, int t, prng_t PRNG) {
	int * a = (int *) calloc(n, sizeof (int));
	for (int i = 0; i < t; ++i) a[i] = 1;
	rand_shuffle(a, n, t, PRNG);
	return a;
}

// mot aléatoire de poids w
f3_t * rand_word(int n, int w, prng_t PRNG) {
	int * perm = randperm(n, PRNG);
	f3_t * e = (f3_t *) malloc(n * sizeof (f3_t));
	for (int i = 0; i < w; ++i) {
		e[perm[i]] = f3_randnonzero(PRNG);
	}
	for (int i = w; i < n; i++) {
		e[perm[i]] = f3_zero();
	}
	return e;
}
