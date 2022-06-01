#ifndef WAVE_H
#define WAVE_H

#include <stdio.h>
#include <stdint.h>

/* Randomness */

typedef struct prng {
	unsigned char buff[64], buff2, buff3, available, available2, available3;
	unsigned long bytecount;
	char * init_str;
} * prng_t;

uint8_t rnd8(prng_t PRNG);
uint8_t rnd_bit(prng_t PRNG);
uint8_t rnd_trit(prng_t PRNG);
uint16_t rnd16(prng_t PRNG);
uint16_t rnd_short(int n, prng_t PRNG);
prng_t prng_init(char * type, char * params_id, unsigned long seed);
void prng_clear(prng_t PRNG);

void rand_shuffle(int * a, int n, int t, prng_t PRNG);
int * identity_perm(int n);
int * randperm(int n, prng_t PRNG);
int * rand_bincw(int n, int t, prng_t PRNG);

/* Wave parameters, including features for precomputed sampling data */

#define RND_MAX_LEN 16
typedef struct {
	unsigned char size;
	unsigned char value[RND_MAX_LEN];
	prng_t PRNG;
} random_stream_t;

#define DISTRIB_PREC 4
// represents the floating point number sum(value[i] / 256^(i+1+offset)) >= 0
// if offset >= 0 the number is < 1
typedef struct {
	char offset;
	unsigned char value[DISTRIB_PREC];
} fp_t;

// (cumulated) distribution of size possible values value[i], 0<=i<size
// x[i]-x[i-1] is the probability to get value[i]
// x[size-1] = 1
typedef struct {
	uint16_t size;
	uint16_t * value;
	fp_t * x;
} Cdist_t;

// array used for rejection
typedef struct {
	int min, max, prec;
	unsigned char ** proba;
} farray_t;

int pickV(prng_t PRNG);
int acceptV(int t, prng_t PRNG);
int pickU(int t, prng_t PRNG);
int acceptU(int t, int l, prng_t PRNG);

void init_params(char * params_id);
void init_precomp(char * params_id);
void init(char * params_id);
void cleanup();

/* Wave parameters END */

/* Ternary arithmetic, including vectors and matrices */

typedef char f3_t;

#define f3_zero() (0)
#define f3_one() (1)
#define f3_two() (-1)
#define f3_iszero(a) ((a) == f3_zero())
#define f3_isnonzero(a) ((a) != f3_zero())
#define f3_isone(a) ((a) == f3_one())
#define f3_istwo(a) ((a) == f3_two())
#define f3_neg(a) (-(a))
#define f3_mul(a, b) ((a) * (b))

f3_t f3_add(f3_t a, f3_t b);

typedef unsigned long long word_t;
// make sure the popcount builtin function matches with word_t
#define popcount __builtin_popcountll

typedef struct {
	word_t * x0, * x1;
	int length, alloc;
} vf3_t;

typedef struct {
	vf3_t * row;
	int rowdim, coldim;
} mf3_t;

#define WORD_LENGTH (8 * sizeof (word_t))

f3_t f3_randnonzero(prng_t PRNG);
f3_t f3_rand(prng_t PRNG);
f3_t * f3_randarray(int n, prng_t PRNG);
f3_t ** f3_randarray2(int k, int n, prng_t PRNG);
int f3_array_write(f3_t * v, int length, FILE * stream);
int f3_array_read(f3_t * v, int length, FILE * stream);
int f3_array_weight(f3_t * v, int length);
char * f3_array_tostring(f3_t * v, int length);

f3_t vf3_coeff(vf3_t x, int j);
int vf3_coeff_iszero(vf3_t x, int j);
int vf3_coeff_isone(vf3_t x, int j);
int vf3_coeff_istwo(vf3_t x, int j);
void vf3_setcoeff(vf3_t x, int j, f3_t a);

vf3_t vf3_new(int length);
void vf3_free(vf3_t x);
void vf3_setrand(vf3_t x, int length, prng_t PRNG);
vf3_t vf3_rand(int length, prng_t PRNG);
void vf3_set(vf3_t x, vf3_t y);
vf3_t vf3_copy(vf3_t y);
void vf3_set_from_array(vf3_t x, f3_t * v, int length);
void vf3_set_permfrom_array(vf3_t x, f3_t * v, int * perm, int length);
vf3_t vf3_from_array(f3_t * a, int length);
void vf3_zero(vf3_t x);
void vf3_toarray(vf3_t x, f3_t * a, int length);
f3_t * vf3_array(vf3_t x, int length);
vf3_t vf3_concat(vf3_t x, int xlen, vf3_t y, int ylen);
void vf3_print(vf3_t x);
// z <- x * y, componentwise
void vf3_mul(vf3_t z, vf3_t x, vf3_t y);
// add all coefficients of x
f3_t vf3_sum(vf3_t x);
// a <- <x,y>
f3_t vf3_dotprod(vf3_t x, vf3_t y);
// z <- x + y
void vf3_add(vf3_t z, vf3_t x, vf3_t y);
// z <- x - y
void vf3_subtract(vf3_t z, vf3_t x, vf3_t y);
// x <- x + y
void vf3_addto(vf3_t x, vf3_t y);
// x <- x - y
void vf3_subtractfrom(vf3_t x, vf3_t y);

f3_t mf3_coeff(mf3_t M, int i, int j);
int mf3_coeff_iszero(mf3_t M, int i, int j);
int mf3_coeff_isone(mf3_t M, int i, int j);
int mf3_coeff_istwo(mf3_t M, int i, int j);
void mf3_setcoeff(mf3_t M, int i, int j, f3_t a);

mf3_t mf3_new(int rowdim, int coldim);
void mf3_free(mf3_t M);
mf3_t mf3_rand(int rowdim, int coldim, prng_t PRNG);
mf3_t mf3_copy(mf3_t P);
void mf3_set(mf3_t P, mf3_t M);
void mf3_zero(mf3_t M);
void mf3_print(mf3_t M);
int mf3_write(mf3_t * M, FILE * stream);
int mf3_read(mf3_t * M, FILE * stream);
mf3_t mf3_diag(f3_t * a, int length);
mf3_t mf3_perm(int * perm, int length);
mf3_t mf3_from_array2(f3_t ** A, int rowdim, int coldim);
mf3_t mf3_permfrom_array2(f3_t ** A, int * perm, int rowdim, int coldim);
mf3_t mf3_augment(mf3_t H, f3_t * s);
int mf3_gauss_elim_single(mf3_t M, int r, int j);
int mf3_gauss_elim(mf3_t M, int * support);
int mf3_gauss_elim_split(mf3_t M, int * support, int t, int * mask);
mf3_t mf3_mul(mf3_t P, mf3_t Q);
mf3_t mf3_stack(mf3_t P, mf3_t Q);
mf3_t mf3_concat(mf3_t P, mf3_t Q);
mf3_t mf3_dual(mf3_t G, int * support);
mf3_t mf3_transpose(mf3_t P);
f3_t * mf3_mv_mul(mf3_t M, vf3_t v);
f3_t * mf3_ma_mul(mf3_t M, f3_t * a);
void mf3_neg(mf3_t M);

/* Ternary arithmetic END */

/* Wave key generation */

typedef struct {
	f3_t *a, *b, *c, *d, *x; // x = ab + dc
} coeff_t;

typedef struct {
	mf3_t HU, HV, Sinv;
	coeff_t coeff;
	int * perm;
} wave_sk_t;

typedef mf3_t wave_pk_t;

void keygen(wave_sk_t * sk, wave_pk_t * pk, prng_t PRNG);
void phiUV(f3_t * e, f3_t * eU, f3_t * eV, wave_sk_t sk);
void wave_sk_clear(wave_sk_t sk);
void wave_pk_clear(wave_pk_t pk);
int fwrite_publickey(wave_pk_t * pk, FILE * stream);
int fread_publickey(wave_pk_t * pk, FILE * stream);
int fwrite_secretkey(wave_sk_t * sk, FILE * stream);
int fread_secretkey(wave_sk_t * sk, FILE * stream);

/* Wave keys END */

/* Wave primitives */

int nz_pairs(f3_t * v, int length);
int * supp(f3_t * v, int length, int * weight);

f3_t * sign_wave(f3_t * s, wave_sk_t sk, prng_t PRNG);
int verif(f3_t * s, f3_t * e, wave_pk_t pk);

/* GLOBAL VARIABLES: system parameters */
int n, n2, kU, kV, k, w, d;

/* GLOBAL VARIABLES: precomputed data */
farray_t rV, * rU;
Cdist_t CLPV, * CLPU;

/* GLOBAL VARIABLES: statistics, others, ... */
unsigned long rejU, rejV, gaussU, gaussV;

int verbose;

/* GLOBAL VARIABLES END */

/* structure for compact parameters description in "params.h" */
typedef struct {
	char * identifier;
	int n, w, kU, kV, d;
} wave_params_t;

#endif
