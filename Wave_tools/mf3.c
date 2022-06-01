/**
 * \file mf3.c : files from Wave implementation to use ternary bitsliced arithmetic
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "wave.h"

f3_t mf3_coeff(mf3_t M, int i, int j)
{
	return vf3_coeff(M.row[i], j);
}

int mf3_coeff_iszero(mf3_t M, int i, int j)
{
	return vf3_coeff_iszero(M.row[i], j);
}

int mf3_coeff_isone(mf3_t M, int i, int j)
{
	return vf3_coeff_isone(M.row[i], j);
}

int mf3_coeff_istwo(mf3_t M, int i, int j)
{
	return vf3_coeff_istwo(M.row[i], j);
}

void mf3_setcoeff(mf3_t M, int i, int j, f3_t a)
{
	vf3_setcoeff(M.row[i], j, a);
}

mf3_t mf3_new(int rowdim, int coldim)
{
	int i;
	mf3_t M;
	M.rowdim = rowdim;
	M.coldim = coldim;
	M.row = (vf3_t *)malloc(M.rowdim * sizeof(vf3_t));
	for (i = 0; i < M.rowdim; ++i)
	{
		M.row[i] = vf3_new(M.coldim);
	}
	return M;
}

void mf3_free(mf3_t M)
{
	int i;
	for (i = 0; i < M.rowdim; ++i)
	{
		vf3_free(M.row[i]);
	}
	free(M.row);
}

mf3_t mf3_copy(mf3_t P)
{
	int i;
	mf3_t M;
	M.rowdim = P.rowdim;
	M.coldim = P.coldim;
	M.row = (vf3_t *)malloc(M.rowdim * sizeof(vf3_t));
	for (i = 0; i < M.rowdim; ++i)
	{
		M.row[i] = vf3_copy(P.row[i]);
	}
	return M;
}

void mf3_set(mf3_t M, mf3_t P)
{
	int i;
	for (i = 0; i < M.rowdim; ++i)
	{
		vf3_set(M.row[i], P.row[i]);
	}
}

void mf3_zero(mf3_t M)
{
	int i;
	for (i = 0; i < M.rowdim; ++i)
	{
		vf3_zero(M.row[i]);
	}
}

void mf3_print(mf3_t M)
{
	int i;
	for (i = 0; i < M.rowdim; ++i)
	{
		vf3_print(M.row[i]);
	}
	printf("\n");
}

int mf3_write(mf3_t *M, FILE *stream)
{
	int i, j;
	f3_t *v = malloc(M->coldim * M->rowdim);
	for (i = 0; i < M->rowdim; ++i)
	{
		vf3_toarray(M->row[i], v + (i * M->coldim), M->coldim);
	}
	j = f3_array_write(v, M->coldim * M->rowdim, stream);
	free(v);
	return j;
}

int mf3_read(mf3_t *M, FILE *stream)
{
	int j = 0;
	f3_t *v = malloc(M->coldim * M->rowdim);
	j = f3_array_read(v, M->coldim * M->rowdim, stream);
	for (int i = 0; i < M->rowdim; ++i)
	{
		vf3_set_from_array(M->row[i], v + (i * M->coldim), M->coldim);
	}
	free(v);
	return j;
}

mf3_t mf3_from_array2(f3_t **A, int rowdim, int coldim)
{
	int i;
	mf3_t M = mf3_new(rowdim, coldim);
	for (i = 0; i < rowdim; ++i)
	{
		vf3_set_from_array(M.row[i], A[i], coldim);
	}
	return M;
}

// Allocate a new matrix Hp equal to H with an additional column equal
// to s.
mf3_t mf3_augment(mf3_t H, f3_t *s)
{
	mf3_t Hp = mf3_new(H.rowdim, H.coldim + 1);
	int i = H.coldim / WORD_LENGTH;
	word_t z = ((word_t)1) << (H.coldim % WORD_LENGTH);
	for (int l = 0; l < H.rowdim; ++l)
	{
		vf3_set(Hp.row[l], H.row[l]);
		if (f3_isone(s[l]))
			Hp.row[l].x0[i] ^= z;
		else if (f3_istwo(s[l]))
			Hp.row[l].x1[i] ^= z;
	}
	// pretend that the last column isn't there
	Hp.coldim = H.coldim; // instead of H.coldim + 1
	return Hp;
}

// single pivot at column j for rank r
int mf3_gauss_elim_single(mf3_t M, int r, int j)
{
	int i;
	f3_t a;
	for (i = r; i < M.rowdim; ++i)
	{
		a = mf3_coeff(M, i, j);
		if (f3_isnonzero(a))
		{
			if (f3_istwo(a))
			{ // normalize row
				// exchanging x0 and x1 multiplies the row by -1
				word_t *temp = M.row[i].x0;
				M.row[i].x0 = M.row[i].x1;
				M.row[i].x1 = temp;
			}
			vf3_t temp = M.row[i];
			M.row[i] = M.row[r];
			M.row[r] = temp;
			for (i = 0; i < M.rowdim; ++i)
			{
				if (i != r)
				{
					a = mf3_coeff(M, i, j);
					if (f3_isone(a))
						vf3_subtractfrom(M.row[i], M.row[r]);
					else if (f3_istwo(a))
						vf3_addto(M.row[i], M.row[r]);
				}
			}
			return 1;
		}
	}
	return 0;
}

/***
		support[] is a permutation of {0,..,n-1}.

		The function returns the total number of pivots tried, successful
		or not.

		After the call support[] will have the r=rank(M) ordered
		successful pivots positions first. The last n-r entries of
		support[] contain the n-r other positions.
***/
int mf3_gauss_elim(mf3_t M, int *support)
{
	int j, k, n, r, d;
	n = M.coldim;
	k = M.rowdim;
	int *pivots = (int *)malloc(k * sizeof(int));
	int *nonpivots = (int *)malloc(n * sizeof(int));
	for (r = 0, d = 0, j = 0; (j < n) && (r < k); ++j)
	{
		if (mf3_gauss_elim_single(M, r, support[j]) == 0)
		{
			nonpivots[d] = support[j];
			d++;
		}
		else
		{
			pivots[r] = support[j];
			r++;
		}
	}
	memcpy(support, pivots, r * sizeof(int));
	memcpy(support + r, nonpivots, d * sizeof(int));
	free(pivots);
	free(nonpivots);
	return j;
}

/***
		support[] is a permutation of {0,..,n-1}.

		If mask[i] != 0, the i-th pivot must be in the first t positions
		of support[] else it must be in the last (n-t) positions of
		support[] the pivots are chosen in the order given by support[]
		and mask[].

		The function returns the total number of pivots tried, successful
		or not. Note that the function exits when the elimination is
		complete or when one of the constraints above is not met.

		After the call support[] will have the r ordered successful pivots
		positions first. The last n-r entries of support[] contain the n-r
		other positions.
***/
int mf3_gauss_elim_split(mf3_t M, int *support, int t, int *mask)
{
	int j0, j1, k, n, r, d;
	n = M.coldim;
	k = M.rowdim;
	int *pivots = (int *)malloc(k * sizeof(int));
	int *nonpivots = (int *)malloc(n * sizeof(int));
	for (r = 0, d = 0, j0 = 0, j1 = 0; r < k; ++r)
	{
		if (mask[r])
		{
			while ((j0 < t) && (mf3_gauss_elim_single(M, r, support[j0]) == 0))
			{
				nonpivots[d] = support[j0];
				++d;
				++j0;
			}
			if (j0 >= t)
				break;
			pivots[r] = support[j0];
			++j0;
		}
		else
		{
			while ((j1 + t < n) && (mf3_gauss_elim_single(M, r, support[j1 + t]) == 0))
			{
				nonpivots[d] = support[j1 + t];
				++d;
				++j1;
			}
			if (j1 + t >= n)
				break;
			pivots[r] = support[j1 + t];
			++j1;
		}
		// invariant: j0 + j1 = d + r
	}
	memcpy(nonpivots + d, support + j0, (t - j0) * sizeof(int));
	d += t - j0;
	memcpy(nonpivots + d, support + (j1 + t), (n - (j1 + t)) * sizeof(int));
	d += n - (j1 + t); // should be n - r
	memcpy(support, pivots, r * sizeof(int));
	memcpy(support + r, nonpivots, d * sizeof(int));
	free(pivots);
	free(nonpivots);
	return j0 + j1;
}

f3_t *mf3_mv_mul(mf3_t M, vf3_t v)
{
	int i;
	f3_t *a = (f3_t *)malloc(M.rowdim * sizeof(f3_t));
	for (i = 0; i < M.rowdim; ++i)
	{
		a[i] = vf3_dotprod(M.row[i], v);
	}
	return a;
}

f3_t *mf3_ma_mul(mf3_t M, f3_t *a)
{
	vf3_t v = vf3_from_array(a, M.row[0].length);
	f3_t *b = mf3_mv_mul(M, v);
	vf3_free(v);
	return b;
}
