/**
 * \file keys.c : files from Wave implementation to use Wave secret keys. 
 */

#include <stdlib.h>
#include <stdio.h>
#include "wave.h"

// e_u, e_v --> e = (e_l, e_r)
void phiUV(f3_t *e, f3_t *eU, f3_t *eV, wave_sk_t sk)
{
	int i, length = sk.HU.coldim;
	for (i = 0; i < length; ++i)
	{
		e[i] = f3_add(f3_mul(sk.coeff.a[i], eU[i]), f3_mul(sk.coeff.b[i], eV[i]));
		e[i + length] = f3_add(f3_mul(sk.coeff.c[i], eU[i]), f3_mul(sk.coeff.d[i], eV[i]));
	}
}

f3_t **full_parity_check_matrix(f3_t **HU, f3_t **HV, coeff_t coeff)
{
	int i, j;
	f3_t **Hsec;
	Hsec = (f3_t **)malloc((n - k) * sizeof(f3_t *));
	for (i = 0; i < n2 - kU; ++i)
	{
		Hsec[i] = (f3_t *)malloc(n * sizeof(f3_t));
		for (j = 0; j < n2; ++j)
		{
			Hsec[i][j] = f3_mul(HU[i][j], coeff.d[j]);
			Hsec[i][j + n2] = f3_neg(f3_mul(HU[i][j], coeff.b[j]));
		}
	}
	for (i = 0; i < n2 - kV; ++i)
	{
		Hsec[i + n2 - kU] = (f3_t *)malloc(n * sizeof(f3_t));
		for (j = 0; j < n2; ++j)
		{
			Hsec[i + n2 - kU][j] = f3_neg(f3_mul(HV[i][j], coeff.c[j]));
			Hsec[i + n2 - kU][j + n2] = f3_mul(HV[i][j], coeff.a[j]);
		}
	}

	return Hsec;
}

void wave_sk_clear(wave_sk_t sk)
{
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

void wave_pk_clear(wave_pk_t pk)
{
	mf3_free(pk);
}

int fwrite_publickey(wave_pk_t *pk, FILE *stream)
{
	return mf3_write(pk, stream);
}

int fread_publickey(wave_pk_t *pk, FILE *stream)
{
	int j = 0;
	*pk = mf3_new(n - k, k);
	j += mf3_read(pk, stream);
	return j;
}

int int_array_write(int *a, int length, FILE *stream)
{
	return sizeof(int) * fwrite(a, sizeof(int), length, stream);
}

int int_array_read(int *a, int length, FILE *stream)
{
	return sizeof(int) * fread(a, sizeof(int), length, stream);
}

int fwrite_secretkey(wave_sk_t *sk, FILE *stream)
{
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

int fread_secretkey(wave_sk_t *sk, FILE *stream)
{
	int j = 0;

	sk->HU = mf3_new(n2 - kU, n2);
	sk->HV = mf3_new(n2 - kV, n2);
	sk->Sinv = mf3_new(n - k, n - k);
	sk->coeff.a = (f3_t *)malloc(n2 * sizeof(f3_t));
	sk->coeff.b = (f3_t *)malloc(n2 * sizeof(f3_t));
	sk->coeff.c = (f3_t *)malloc(n2 * sizeof(f3_t));
	sk->coeff.d = (f3_t *)malloc(n2 * sizeof(f3_t));
	sk->coeff.x = (f3_t *)malloc(n2 * sizeof(f3_t));
	sk->perm = (int *)malloc(n * sizeof(int));

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

void print_publickey(wave_pk_t *pk)
{
	printf("Public key:\n");
	mf3_print(*pk);
}

void print_secretkey(wave_sk_t *sk, int verbose)
{
	printf("Secret key:\n");
	if (verbose > 0)
	{
		printf("HU:\n");
		mf3_print(sk->HU);
		printf("HV:\n");
		mf3_print(sk->HV);
		printf("Sinv:\n");
		mf3_print(sk->Sinv);
	}
	printf("coeff: a : %d, b : %d, c : %d, d : %d\n", (int)*sk->coeff.a, (int)*sk->coeff.b, (int)*sk->coeff.c, (int)*sk->coeff.d);
	printf("perm:\n");
	int i;
	for (i = 0; i < ((verbose > 0) ? n : 10); ++i)
	{
		printf("pi(%d) = %d ", i, sk->perm[i]);
	}
	printf("\n");
}