/**
 * \file mf3.c : files from Wave implementation to use ternary bitsliced arithmetic
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "wave.h"

/*** ternary field elements ***/
f3_t f3_add(f3_t a, f3_t b)
{
	static f3_t f3add[5] = {1, -1, 0, 1, -1};
	return f3add[2 + a + b];
}

void base3(unsigned char c, f3_t *d, int size)
{
	for (int i = 0; i < size; ++i)
	{
		d[size - 1 - i] = (c % 3) - 1; // d[i] in {-1,0,1} <= {0,1,2}
		c /= 3;
	}
}

unsigned char compact3(f3_t *d, int size)
{
	unsigned char c = 0;
	for (int i = 0; i < size; ++i)
	{
		c *= 3;
		c += (d[i] + 1); // d[i] in {-1,0,1} => {0,1,2}
	}
	return c;
}

int f3_array_write(f3_t *v, int length, FILE *stream)
{
	int i;
	for (i = 0; i < length - 5; i += 5)
	{
		putc(compact3(v + i, 5), stream);
	}
	if (i < length)
	{
		putc(compact3(v + i, length - i), stream);
	}
	return (length - 1) / 5 + 1;
}

int f3_array_read(f3_t *v, int length, FILE *stream)
{
	int i, j, c;
	for (i = 0, j = 0; i < length - 5; i += 5)
	{
		c = getc(stream);
		if (c == EOF)
			return j;
		base3(c, v + i, 5);
		j += 5;
	}
	if (i < length)
	{
		c = getc(stream);
		if (c == EOF)
			return j;
		base3(c, v + i, length - i);
		j += (length - i);
	}
	return j;
}

// Hamming weight of v
int f3_array_weight(f3_t *v, int length)
{
	int i, j;
	for (i = 0, j = 0; i < length; ++i)
		if (f3_isnonzero(v[i]))
			++j;
	return j;
}

void print_f3(f3_t a)
{
	printf("%d", a);
}

/***
		ternary field arrays (bitsliced)

		the element a in {0, 1, -1} is represented by a pair of bits (a0, a1)
		in {(0, 0), (1, 0), (0, 1)} such that a = a0 - a1

		the pair of binary words (x0, x1) will represent a ternary array,
		the i-th element of this array will be x0i - x1i where x0i and
		x1i are the i-th bit of x0 and x1 respectively
 ***/
vf3_t vf3_new(int length)
{
	vf3_t x;
	x.length = length;
	x.alloc = 1 + (length - 1) / WORD_LENGTH;
	x.x0 = (word_t *)calloc(x.alloc, sizeof(word_t));
	x.x1 = (word_t *)calloc(x.alloc, sizeof(word_t));
	return x;
}

void vf3_free(vf3_t x)
{
	free(x.x0);
	free(x.x1);
}

f3_t vf3_coeff(vf3_t x, int j)
{
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return ((x.x0[i] >> j) & 1) - ((x.x1[i] >> j) & 1);
}

int vf3_coeff_iszero(vf3_t x, int j)
{
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (((x.x0[i] | x.x1[i]) >> j) & 1) == 0;
}

int vf3_coeff_isone(vf3_t x, int j)
{
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (x.x0[i] >> j) & 1;
}

int vf3_coeff_istwo(vf3_t x, int j)
{
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (x.x1[i] >> j) & 1;
}

void vf3_setcoeff(vf3_t x, int j, f3_t a)
{
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	word_t z = ((word_t)1) << j;
	if (f3_isone(a))
	{
		x.x0[i] |= z;
		x.x1[i] &= ~z;
	}
	else if (f3_istwo(a))
	{
		x.x0[i] &= ~z;
		x.x1[i] |= z;
	}
	else
	{
		x.x0[i] &= ~z;
		x.x1[i] &= ~z;
	}
}

// assumes x.alloc >= y.alloc
void vf3_set(vf3_t x, vf3_t y)
{
	memcpy(x.x0, y.x0, y.alloc * sizeof(word_t));
	memcpy(x.x1, y.x1, y.alloc * sizeof(word_t));
}

vf3_t vf3_copy(vf3_t y)
{
	vf3_t x = vf3_new(y.length);
	vf3_set(x, y);
	return x;
}

void vf3_set_from_array(vf3_t x, f3_t *v, int length)
{
	int i, j;
	word_t z;

	vf3_zero(x);
	for (i = 0, j = 0, z = 1; j < length; ++j, z <<= 1)
	{
		if (z == 0)
		{
			z = 1;
			++i;
		}
		if (f3_isone(v[j]))
			x.x0[i] ^= z;
		else if (f3_istwo(v[j]))
			x.x1[i] ^= z;
	}
}

vf3_t vf3_from_array(f3_t *a, int length)
{
	vf3_t v = vf3_new(length);
	vf3_set_from_array(v, a, length);
	return v;
}

void vf3_zero(vf3_t x)
{
	memset(x.x0, 0, x.alloc * sizeof(word_t));
	memset(x.x1, 0, x.alloc * sizeof(word_t));
}

void vf3_toarray(vf3_t x, f3_t *a, int length)
{
	int i, j;
	word_t z;

	for (i = 0, j = 0, z = 1; j < length; ++j, z <<= 1)
	{
		if (z == 0)
		{
			z = 1;
			++i;
		}
		if (x.x0[i] & z)
			a[j] = f3_one();
		else if (x.x1[i] & z)
			a[j] = f3_two();
		else
			a[j] = f3_zero();
	}
}

f3_t *vf3_array(vf3_t x, int length)
{
	f3_t *a = (f3_t *)malloc(length * sizeof(f3_t));

	vf3_toarray(x, a, length);

	return a;
}

void vf3_print(vf3_t x)
{
	int i, j;
	word_t z;

	for (i = 0, j = 0, z = 1; j < x.length; ++j, z <<= 1)
	{
		if (z == 0)
		{
			z = 1;
			++i;
		}
		int a = 0;
		if (x.x0[i] & z)
			a = 1;
		else if (x.x1[i] & z)
			a = 2;
		printf("%d", a);
	}
	printf("\n");
}

/*
	 bitsliced arithmetic for multiplying (componentwise) ternary vectors

	 (a0, a1) in {(0, 0), (1, 0), (0,1)} represent respectively 0, 1,
	 and -1, that is a = a0 - a1

	 r0, r1, x0, x1, y0, y1 are binary and represent r, x, and y
	 to get r = x * y, we want to express r0 and r1 such that
	 r0 - r1 = (x0 - x1) * (y0 - y1)

	 modulo 2 we have
	 r0 = (x0 * y0) + (x1 * y1)
	 r1 = (x0 * y1) + (x1 * y0)

	 which can be computed with 6 logical operations (see below)
	 open problem: do better or prove that this is optimal
	 (i.e. find other, easier to compute, identities for r0 and r1 or
	 find a better way to evaluate the identities above)
*/
void vf3_mul_aux(word_t *r0, word_t *r1, word_t x0, word_t x1, word_t y0, word_t y1)
{
	word_t a, b, c, d;
	a = x0 & y0;
	b = x1 & y1;
	c = x0 & y1;
	d = x1 & y0;
	*r0 = a ^ b;
	*r1 = c ^ d;
}

// z <- x * y
// componentwise product
void vf3_mul(vf3_t z, vf3_t x, vf3_t y)
{
	int i;
	for (i = 0; i < x.alloc; ++i)
	{
		vf3_mul_aux(z.x0 + i, z.x1 + i, x.x0[i], x.x1[i], y.x0[i], y.x1[i]);
	}
}

int binary_weight(word_t *x, int n)
{
	int i, w = 0;
	for (i = 0; i < n; ++i)
		w += popcount(x[i]);
	return w;
}

f3_t vf3_sum(vf3_t x)
{
	// f3_t = {0,1,-1}
	int j = (binary_weight(x.x0, x.alloc) + 2 * binary_weight(x.x1, x.alloc)) % 3;
	return (j == 2) ? (-1) : j;
}

f3_t vf3_dotprod(vf3_t x, vf3_t y)
{
	vf3_t z = vf3_new(x.length);
	vf3_mul(z, x, y);
	f3_t a = vf3_sum(z);
	vf3_free(z);
	return a;
}

/*
	 bitsliced arithmetic for adding ternary vectors

	 (a0, a1) in {(0, 0), (1, 0), (0,1)} represents respectively 0, 1,
	 and -1, that is a = a0 - a1

	 (r0, r1), (x0, x1), (y0, y1) are binary and represent r, x, y
	 to get r = x + y, we want to express r0 and r1 so that we always have
	 r0 - r1 = (x0 - x1) + (y0 - y1)

	 modulo 2 we have
	 r0 = (x0 + y0) * (x1 + y1) + (x1 * y1) + (x0 + y0)
	 r1 = (x0 + y0) * (x1 + y1) + (x0 * y0) + (x1 + y1)

	 which can be computed with 9 logical operations (see below)
	 open problem: do better or prove that this is optimal
	 (i.e. find other, easier to compute, identities for r0 and r1 or
	 find a better way to evaluate the identities above)
*/
void vf3_add_aux(word_t *r0, word_t *r1, word_t x0, word_t x1, word_t y0, word_t y1)
{
	word_t a, b, c, d, e;
	a = x0 ^ y0;
	b = x1 ^ y1;
	c = x0 & y0;
	d = x1 & y1;
	e = a & b;
	*r0 = e ^ d ^ a;
	*r1 = e ^ b ^ c;
}

// z <- x + y
void vf3_add(vf3_t z, vf3_t x, vf3_t y)
{
	int i;
	for (i = 0; i < x.alloc; ++i)
	{
		vf3_add_aux(z.x0 + i, z.x1 + i, x.x0[i], x.x1[i], y.x0[i], y.x1[i]);
	}
}

// z <- x - y
void vf3_subtract(vf3_t z, vf3_t x, vf3_t y)
{
	int i;
	for (i = 0; i < x.alloc; ++i)
	{
		vf3_add_aux(z.x0 + i, z.x1 + i, x.x0[i], x.x1[i], y.x1[i], y.x0[i]);
	}
}

// x <- x + y
void vf3_addto(vf3_t x, vf3_t y)
{
	int i;
	for (i = 0; i < x.alloc; ++i)
	{
		vf3_add_aux(x.x0 + i, x.x1 + i, x.x0[i], x.x1[i], y.x0[i], y.x1[i]);
	}
}

// x <- x - y
void vf3_subtractfrom(vf3_t x, vf3_t y)
{
	int i;
	for (i = 0; i < x.alloc; ++i)
	{
		vf3_add_aux(x.x0 + i, x.x1 + i, x.x0[i], x.x1[i], y.x1[i], y.x0[i]);
	}
}
