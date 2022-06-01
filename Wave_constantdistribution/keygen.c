#include <stdlib.h>
#include <stdio.h>
#include "wave.h"

coeff_t coeffUV(int n, prng_t PRNG) {
	int i;
	coeff_t coeff;
	coeff.a = (f3_t *) malloc(n * sizeof (f3_t));
	coeff.b = (f3_t *) malloc(n * sizeof (f3_t));
	coeff.c = (f3_t *) malloc(n * sizeof (f3_t));
	coeff.d = (f3_t *) malloc(n * sizeof (f3_t));
	coeff.x = (f3_t *) malloc(n * sizeof (f3_t));

	for (i = 0; i < n; ++i) {
		coeff.a[i] = f3_randnonzero(PRNG);
		coeff.c[i] = f3_randnonzero(PRNG);
		coeff.b[i] = f3_rand(PRNG);
		coeff.d[i] = f3_mul(coeff.a[i], f3_add(f3_one(), f3_mul(coeff.b[i], coeff.c[i])));
		coeff.x[i] = f3_add(f3_mul(coeff.a[i], coeff.b[i]), f3_mul(coeff.c[i], coeff.d[i]));
	}

	return coeff;
}

// e_u, e_v --> e = (e_l, e_r)
void phiUV(f3_t * e, f3_t * eU, f3_t * eV, wave_sk_t sk) {
	int i, length = sk.HU.coldim;
	for (i = 0; i < length; ++i) {
		e[i] = f3_add(f3_mul(sk.coeff.a[i], eU[i]), f3_mul(sk.coeff.b[i], eV[i]));
		e[i + length] = f3_add(f3_mul(sk.coeff.c[i], eU[i]), f3_mul(sk.coeff.d[i], eV[i]));
	}
}

f3_t ** full_parity_check_matrix(f3_t ** HU, f3_t ** HV, coeff_t coeff) {
	int i, j;
	f3_t ** Hsec;
	Hsec = (f3_t **) malloc((n - k) * sizeof (f3_t *));
	for (i = 0; i < n2 - kU; ++i) {
		Hsec[i] = (f3_t *) malloc(n * sizeof (f3_t));
		for (j = 0; j < n2; ++j) {
			Hsec[i][j] = f3_mul(HU[i][j], coeff.d[j]);
			Hsec[i][j + n2] = f3_neg(f3_mul(HU[i][j], coeff.b[j]));
		}
	}
	for (i = 0; i < n2 - kV; ++i) {
		Hsec[i + n2 - kU] = (f3_t *) malloc(n * sizeof (f3_t));
		for (j = 0; j < n2; ++j) {
			Hsec[i + n2 - kU][j] = f3_neg(f3_mul(HV[i][j], coeff.c[j]));
			Hsec[i + n2 - kU][j + n2] = f3_mul(HV[i][j], coeff.a[j]);
		}			
	}

	return Hsec;
}

void keygen(wave_sk_t * sk, wave_pk_t * pk, prng_t PRNG) {
	int i, j;
	f3_t ** Hsec;

	// random permutation from PRNG
	int * perm = randperm(n, PRNG);

	sk->perm = perm;

	// random matrices and coefficients from PRNG
	f3_t ** HU = f3_randarray2(n2 - kU, n2, PRNG);
	f3_t ** HV = f3_randarray2(n2 - kV, n2, PRNG);
	coeff_t coeff = coeffUV(n2, PRNG);

	sk->HU = mf3_from_array2(HU, n2 - kU, n2);
	sk->HV = mf3_from_array2(HV, n2 - kV, n2);
	sk->coeff = coeff;

	// public key
	// first, build the full (secret) parity check matrix H
	Hsec = full_parity_check_matrix(HU, HV, coeff);
	mf3_t H = mf3_from_array2(Hsec, n - k, n);
	for (i = 0; i < n2 - kU; ++i) {
		free(HU[i]);
	}
	free(HU);
	for (i = 0; i < n2 - kV; ++i) {
		free(HV[i]);
	}
	free(HV);

	// second, Gaussian elimination, choosing the pivot position in the
	// order given by perm
	mf3_gauss_elim(H, perm); // d positions sont modifié et avec une très grande probabilité d < gap. 
	// perm is possibly modified such that its first (n-k) entries give
	// the actual pivot positions

	// finally, extract the public key from H (the redondancy part)
	// according to the secret permutation perm, the last k entries of
	// perm contain the non pivot positions
	mf3_t R = mf3_new(n - k, k);
	f3_t * v = (f3_t *) malloc(k * sizeof (f3_t));
	for (i = 0; i < n - k; ++i) {
		for (j = 0; j < k; ++j) {
			v[j] = mf3_coeff(H, i, perm[n - k + j]);
		}
		vf3_set_from_array(R.row[i], v, k);
	}
	free(v);
	mf3_free(H);

	// additional step, we store in the secret key the (n-k)x(n-k)
	// matrix given by the first (n-k) columns (in the order given by
	// perm) of Hsec
	mf3_t Sinv = mf3_new(n - k, n - k);
	v = (f3_t *) malloc((n - k) * sizeof (f3_t));
	for (i = 0; i < n - k; ++i) {
		for (j = 0; j < n - k; ++j) {
			v[j] = Hsec[i][perm[j]];
		}
		vf3_set_from_array(Sinv.row[i], v, n - k);
	}
	for (i = 0; i < n - k; ++i) {
		free(Hsec[i]);
	}
	free(Hsec);
	free(v);
	sk->Sinv = Sinv;
	*pk = R;
}

void wave_sk_clear(wave_sk_t sk) {
	free(sk.coeff.a);
	free(sk.coeff.b);
	free(sk.coeff.c);
	free(sk.coeff.d);
	free(sk.coeff.x);
	free(sk.perm);
	mf3_free(sk.HU);
	mf3_free(sk.HV);
	mf3_free(sk.Sinv);
}

void wave_pk_clear(wave_pk_t pk) {
	mf3_free(pk);
}

int fwrite_publickey(wave_pk_t * pk, FILE * stream) {
	return mf3_write(pk, stream);
}

int fread_publickey(wave_pk_t * pk, FILE * stream) {
	int j = 0;
	*pk = mf3_new(n - k, k);
	j += mf3_read(pk, stream);
	return j;
}

int int_array_write(int * a, int length, FILE * stream) {
	return sizeof (int) * fwrite(a, sizeof (int), length, stream);
}

int int_array_read(int * a, int length, FILE * stream) {
	return sizeof (int) * fread(a, sizeof (int), length, stream);
}

int fwrite_secretkey(wave_sk_t * sk, FILE * stream) {
	int j = 0;
	j += mf3_write(&(sk->HU), stream);
	j += mf3_write(&(sk->HV), stream);
	j += mf3_write(&(sk->Sinv), stream);
	j += f3_array_write(sk->coeff.a, n2, stream);
	j += f3_array_write(sk->coeff.b, n2, stream);
	j += f3_array_write(sk->coeff.c, n2, stream);
	j += f3_array_write(sk->coeff.d, n2, stream);
	j += f3_array_write(sk->coeff.x, n2, stream);
	j += int_array_write(sk->perm, n, stream);
	return j;
}

int fread_secretkey(wave_sk_t * sk, FILE * stream) {
	int j = 0;

	sk->HU = mf3_new(n2 - kU, n2);
	sk->HV = mf3_new(n2 - kV, n2);
	sk->Sinv = mf3_new(n - k, n - k);
	sk->coeff.a = (f3_t *) malloc(n2 * sizeof (f3_t));
	sk->coeff.b = (f3_t *) malloc(n2 * sizeof (f3_t));
	sk->coeff.c = (f3_t *) malloc(n2 * sizeof (f3_t));
	sk->coeff.d = (f3_t *) malloc(n2 * sizeof (f3_t));
	sk->coeff.x = (f3_t *) malloc(n2 * sizeof (f3_t));
	sk->perm = (int *) malloc(n * sizeof (int));

	j += mf3_read(&(sk->HU), stream);
	j += mf3_read(&(sk->HV), stream);
	j += mf3_read(&(sk->Sinv), stream);
	j += f3_array_read(sk->coeff.a, n2, stream);
	j += f3_array_read(sk->coeff.b, n2, stream);
	j += f3_array_read(sk->coeff.c, n2, stream);
	j += f3_array_read(sk->coeff.d, n2, stream);
	j += f3_array_read(sk->coeff.x, n2, stream);
	j += int_array_read(sk->perm, n, stream);

	return j;
}
