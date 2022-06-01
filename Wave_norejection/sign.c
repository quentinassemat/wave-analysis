#include <stdlib.h>
#include <string.h>
#include "wave.h"

// v a vector of size 2 * length
// number of non-zero pairs (v[i], v[i+length])
int nz_pairs(f3_t * v, int length) {
	int i, l;
	for (i = 0, l = 0; i < length; ++i)
		if (f3_iszero(v[i]) && f3_iszero(v[i + length]))
			++l;
	return length - l;
}

// v is a vector of size length.
// support will contain all integers from 0 to length-1.
// Non-zero postions of v come first in support.
// Side effect: weight is set to the Hamming weight of v.
int * supp(f3_t * v, int length, int * weight) {
	int i, j, l;
	int * support = malloc(length * sizeof (int));
	for (j = 0, i = 0, l = length - 1; j < length; ++j) {
		if (f3_isnonzero(v[j])) {
			support[i] = j;
			++i;
		}
		else {
			support[l] = j;
			--l;
		}
	}
	if (weight != NULL)
		*weight = i;
	return support;
}

f3_t * decodeV(f3_t * sV, mf3_t HV, prng_t PRNG) {
	int i, t;
	f3_t * y;
	int * support = identity_perm(n2);
	mf3_t H = mf3_augment(HV, sV);
	f3_t * eV = (f3_t *) malloc((n2 + 1) * sizeof (f3_t));
	eV[n2] = f3_neg(f3_one()); // -1
	do {
		rand_shuffle(support, n2, n2 - kV + d, PRNG);
		gaussV++;
	} while (mf3_gauss_elim(H, support) > n2 - kV + d);
	/*
		mf3_gauss_elim() returns the number of pivots that were tried
		before the Guassian elimination ended.

		If d is large enough (the fully proven variant of Wave) the
		above condition succeeds with probability overwhelmingly close to 1.

		If d = 0 (the 'no gap' variant) we repeat the Gaussian
		elimination until no pivot fails.
	*/
	while (1) {
		t = pickV(PRNG);
		for (i = 0; i < n2 - kV; ++i)
			eV[support[i]] = f3_zero();
		for (; i < n2 - kV + d; ++i)
			eV[support[i]] = f3_rand(PRNG);
		for (; i < n2 - t; ++i)
			eV[support[i]] = f3_zero();
		for (; i < n2; ++i)
			eV[support[i]] = f3_randnonzero(PRNG);
		y = mf3_ma_mul(H, eV);
		for (i = 0; i < n2 - kV; ++i)
			eV[support[i]] = f3_neg(y[i]);
		free(y);
		if (acceptV(f3_array_weight(eV, n2), PRNG))
			break;
		rejV++;
	}
	free(support);
	mf3_free(H);
	return eV;
}

f3_t * decodeU(f3_t * sU, f3_t * eV, mf3_t HU, wave_sk_t sk, prng_t PRNG) {
	int i, t, k_nz;
	f3_t * y;
	int * suppV = supp(eV, n2, &t);
	/*
		supp(eV, n2, &t) returns an array of all coordinates
		{0, 1 ... n2-1} starting with the non zero positions of eV.
		Side effect: the value of t is set to the weight of eV.
	*/
	int * support = (int *) malloc(n2 * sizeof (int));
	mf3_t H = mf3_augment(HU, sU);
	f3_t * e = (f3_t *) malloc(n2 * 2 * sizeof (f3_t));
	f3_t * eU = (f3_t *) malloc((n2 + 1) * sizeof (f3_t));
	eU[n2] = f3_neg(f3_one()); // -1
	while (1) {
		k_nz = pickU(t, PRNG);
		do {
			rand_shuffle(suppV, t, t - k_nz, PRNG);
			rand_shuffle(suppV + t, n2 - t, n2 - kU + d - (t - k_nz), PRNG);
			memcpy(support, suppV, (t - k_nz) * sizeof (int));
			memcpy(support + (t - k_nz), suppV + t, (n2 - t) * sizeof (int));
			memcpy(support + (n2 - k_nz), suppV + (t - k_nz), k_nz * sizeof (int));
			/*
				Before the Gaussian elimination, support[] has
				1- (t - k_nz) elements of support(e_V) chosen uniformly at
				random in its first n2 - kU + d entries
				2- (n2 - kU + d - (t - k_nz)) elements out of support(e_V)
				chosen uniformly at random in its first n2 - kU + d
				entries
			 */
			gaussU++;
		} while (mf3_gauss_elim(H, support) > n2 - kU + d);
		/*
			After the Gaussian elimination, support[] still verifies the
			above two conditions and in addition
			3- its first (n2 - kU) entries are the pivots positions

			If d is large enough (the fully proven variant of Wave) the
			elimination succeed with probability overwhelmingly close to 1.

			If d = 0 (the 'no gap' variant) we repeat the Gaussian
			elimination until no pivot fails.
		*/
		do {
			for (i = 0; i < n2 - kU; ++i)
				eU[support[i]] = f3_zero();
			for (; i < n2 - kU + d; ++i)
				eU[support[i]] = f3_rand(PRNG);
			for (; i < n2; ++i)
				if (f3_iszero(eV[support[i]]))
					eU[support[i]] = f3_randnonzero(PRNG);
				else
					eU[support[i]] = f3_mul(sk.coeff.x[support[i]], eV[support[i]]);
			y = mf3_ma_mul(H, eU);
			for (i = 0; i < n2 - kU; ++i)
				eU[support[i]] = f3_neg(y[i]);
			free(y);
			phiUV(e, eU, eV, sk);
		} while (f3_array_weight(e, n) != w);
		if (acceptU(t, nz_pairs(e, n2), PRNG))
			break;
		rejU++;
	}
	free(suppV);
	mf3_free(H);
	free(eU);
	free(support);
	return e;
}

f3_t * sign_wave(f3_t * s, wave_sk_t sk, prng_t PRNG) {
	f3_t * sprime = mf3_ma_mul(sk.Sinv, s);
	f3_t * eV = decodeV(sprime + (n2 - kU), sk.HV, PRNG);
	f3_t * eprime = decodeU(sprime, eV, sk.HU, sk, PRNG);
	f3_t * e = (f3_t *) malloc(n * sizeof (f3_t));
	for (int i = 0; i < n; ++i)
		e[i] = eprime[sk.perm[i]];
	free(sprime);
	free(eprime);
	free(eV);
	return e;
}
