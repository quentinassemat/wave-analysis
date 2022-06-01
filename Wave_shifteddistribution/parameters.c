#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "wave.h"
#include "params.h"

int cmp(random_stream_t * r, fp_t x, int i) {
	if (i - x.offset >= DISTRIB_PREC)
		return 0;
	if (r->size <= i) {
		r->value[i] = rnd8(r->PRNG);
		r->size++;
	}
	unsigned char c = (i - x.offset < 0) ? 0 : x.value[i - x.offset];
	if (c == r->value[i])
		return cmp(r, x, i + 1);
	else if (r->value[i] < c)
		return -1;
	else
		return 1;
}

int pickV(prng_t PRNG) {
	int i_min = -1;
	int i_max = CLPV.size - 1;
	random_stream_t r;
	r.size = 0;
	r.PRNG = PRNG;
	while (i_max - i_min > 1) {
		int i = (i_min + i_max) / 2;
		if (cmp(&r, CLPV.x[i], 0) < 0)
			i_max = i;
		else
			i_min = i;
	}
	return CLPV.value[i_max];
}

// int acceptV(int t, prng_t PRNG) {
// 	if ((t < rV.min) || (t > rV.max))
// 		return 0;
// 	for (int i = 0; i < rV.prec / 8; ++i) {
// 		unsigned char r = rnd8(PRNG);
// 		if (r > rV.proba[t - 1][i])
// 			return 1;
// 		else if (r < rV.proba[t - 1][i])
// 			return 0;
// 	}
// 	return 1;
// }

int acceptV(int t, prng_t PRNG) {
	if ((t < rV.min) || (t > rV.max))
		return 0;
	for (int i = 0; i < rV.prec / 8; ++i) {
		unsigned char r = rnd8(PRNG);
		if (r > rV.proba[t + 1][i])
			return 1;
		else if (r < rV.proba[t + 1][i])
			return 0;
	}
	return 1;
}

int pickU(int t, prng_t PRNG) {
	int i_min = -1;
	int i_max = CLPU[t].size - 1;
	random_stream_t r;
	r.size = 0;
	r.PRNG = PRNG;
	while (i_max - i_min > 1) {
		int i = (i_min + i_max) / 2;
		if (cmp(&r, CLPU[t].x[i], 0) < 0)
			i_max = i;
		else
			i_min = i;
	}
	return CLPU[t].value[i_max];
}

// int acceptU(int t, int l, prng_t PRNG) {
// 	if ((l < rU[t].min) || (l > rU[t].max))
// 		return 0;
// 	for (int i = 0; i < rU[t].prec / 8; ++i) {
// 		unsigned char r = rnd8(PRNG);
// 		if (r > rU[t].proba[l - 1][i])
// 			return 1;
// 		else if (r < rU[t].proba[l - 1][i])
// 			return 0;
// 	}
// 	return 1;
// }

int acceptU(int t, int l, prng_t PRNG) {
	if ((l < rU[t].min) || (l > rU[t].max))
		return 0;
	for (int i = 0; i < rU[t].prec / 8; ++i) {
		unsigned char r = rnd8(PRNG);
		if (r > rU[t].proba[l + 1][i])
			return 1;
		else if (r < rU[t].proba[l + 1][i])
			return 0;
	}
	return 1;
}

Cdist_t read_Cdist_raw24(FILE * stream) {
	int i, j, l;
	uint16_t s;
	unsigned char e, c;
	unsigned u;
	Cdist_t d;

	j = 2 * fread(&s, 2, 1, stream);
	d.size = s;
	d.value = malloc(d.size * sizeof (uint16_t));
	d.x = malloc(d.size * sizeof (fp_t));
	//	printf("%d\n", s);
	for (i = 0; i < d.size; ++i) {
		j += 2 * fread(&s, 2, 1, stream);
		d.value[i] = s;
		e = getc(stream);
		c = getc(stream);
		u = c;
		c = getc(stream);
		u = u ^ (c << 8);
		c = getc(stream);
		u = u ^ (c << 16);
		j += 4;
		u <<= 8 - (e % 8);
		d.x[i].offset = (e / 8) - 3;
		for (l = 0; l < 4; ++l) {
			d.x[i].value[l] = (u >> (24 - l * 8)) & 0xff;
		}
	}
	if (j != 6 * d.size + 2) {
		printf("read error %d %d\n", j, 6 * d.size + 2);
	}

	return d;
}

farray_t read_farray_raw(FILE * stream) {
	int i, l;
	farray_t d;

	// assume d.max < 2^16
	d.min = getc(stream);
	d.min += getc(stream) << 8;
	d.max = getc(stream);
	d.max += getc(stream) << 8;
	d.proba = (unsigned char **) calloc(d.max + 1, sizeof (unsigned char *));
	l = getc(stream);
	d.prec = l * 8; // overides default precision

	unsigned char * u = malloc(l * (d.max - d.min + 1));
	i = fread(u, 1, l * (d.max - d.min + 1), stream);
	for (i = d.min; i <= d.max; ++i, u += l) {		
		d.proba[i] = u;
	}
	return d;
}

char * current_params_id = NULL;
char * precomp_id = NULL;

void init_params(char * params_id) {
	if ((current_params_id == NULL)
			|| (strcmp(current_params_id, params_id) != 0)) {
		current_params_id = realloc(current_params_id, strlen(params_id) + 1);
		strcpy(current_params_id, params_id);
	}
	int i = 0;
	while (1) {
		if (strcmp(defined_params[i].identifier, "") == 0) {
			printf("unknown parameters identifier '%s'\n", params_id);
			exit(0);
		}
		if (strcmp(defined_params[i].identifier, params_id) == 0) {
			n = defined_params[i].n;
			w = defined_params[i].w;
			kU = defined_params[i].kU;
			kV = defined_params[i].kV;
			d = defined_params[i].d;
			break;
		}
		++i;
	}
	n2 = n / 2;
	k = kU + kV;
}

void precomp_cleanup() {
	if (precomp_id == NULL)
		return;
	free(rV.proba[rV.min]);
	free(rV.proba);
	free(CLPV.value);
	free(CLPV.x);
	for (int t = rV.min; t <= rV.max; ++t) {
		free(rU[t].proba[rU[t].min]);
		free(rU[t].proba);
		free(CLPU[t].value);
		free(CLPU[t].x);
	}
	free(rU);
	free(CLPU);
	free(precomp_id);
	precomp_id = NULL;
}

char * precomp_filename(char * params_id) {
	char * filename;

#define PRECOMP_FILE_FORMAT "../Data/%s/wave_precomp%s.dat"
	filename = malloc(strlen(PRECOMP_FILE_FORMAT) + 2 * strlen(params_id));
	sprintf(filename, PRECOMP_FILE_FORMAT, params_id, params_id);
	return filename;
}

void init_precomp(char * params_id) {
	if (precomp_id != NULL) {
		if (strcmp(precomp_id, params_id) == 0)
			return;
		else
			precomp_cleanup();
	}

	precomp_id = realloc(precomp_id, strlen(params_id) + 1);
	strcpy(precomp_id, params_id);

	char * filename = precomp_filename(params_id);
	FILE * stream = fopen(filename, "r");
	if (stream == NULL) {
		printf("Cannot open file %s!\n", filename);
		exit(0);
	}

	rV = read_farray_raw(stream);
	CLPV = read_Cdist_raw24(stream);
	rU = (farray_t *) malloc((rV.max + 1) * sizeof (farray_t));
	CLPU = (Cdist_t *) malloc((rV.max + 1) * sizeof (Cdist_t));
	for (int t = rV.min; t <= rV.max; ++t) {
		rU[t] = read_farray_raw(stream);
		CLPU[t] = read_Cdist_raw24(stream);
	}

	free(filename);
	fclose(stream);
}

void init(char * params_id) {
	init_params(params_id);
	init_precomp(params_id);
}

void cleanup() {
	precomp_cleanup();
	if (current_params_id != NULL)
		free(current_params_id);
	current_params_id = NULL;
}
