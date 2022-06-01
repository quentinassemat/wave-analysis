#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cmdline.h"
#include "wave.h"

extern farray_t rV, * rU;
extern Cdist_t CLPV, * CLPU;

double probaf(farray_t d, int i) {
	int j;
	double x;
	for (j = 0, x = 0; j < d.prec / 8; ++j) {
		x += ((double) d.proba[i][j]) * exp2(-8 * (j + 1));
	}
	return x;
}

void show_farray(farray_t d) {
	int i;

	printf("%d\t%d\n", d.min, d.max);
	for (i = d.min; i <= d.max; ++i) {
		printf("%d\t", i);
		for (int j = 0; j < d.prec / 8; ++j)
			printf("%02x ", d.proba[i][j]);
		printf("\t%g\n", probaf(d, i));
	}
}

double probaC(fp_t x) {
	int j;
	double y;
	for (j = 0, y = 0; j < DISTRIB_PREC; ++j) {
		y += ((double) x.value[j]) * exp2(-8*(j + 1 + x.offset));
	}
	return y;
}

double show_Cdist(Cdist_t d) {
	int i, j, l, b;

	b = rV.prec / 8;

	double y = 0, z = 0;
	printf("%d\n", d.size);
	for (i = 0; i < d.size; ++i) {
		printf("%d\t", d.value[i]);
		double x = probaC(d.x[i]);
		z += (x - y) * d.value[i];
		y = x;
		if (x >= 1) {
			for (j = 0; j < b; ++j) {
				printf("00 ");
			}
		}
		else {
			for (j = d.x[i].offset, l = 0; j > 0; --j, ++l) {
				printf("00 ");
			}
			for (j = 0; (j < DISTRIB_PREC) && (l < b); ++j, ++l) {
				printf("%02x ", d.x[i].value[j]);
			}
			for (; l < b; ++l) {
				printf("00 ");
			}
		}
		printf("\t%g\n", x);
	}
	return z;
}

void show_precomp() {
	printf("internal distribution for V, nb of entries: ");
	show_Cdist(CLPV);
	printf("rejection probabilities for V, range of values: ");
	show_farray(rV);
	for (int t = rV.min; t <= rV.max; ++t) {
		printf("internal distribution for U, |e_V| = %d, nb of entries: ", t);
		double z = show_Cdist(CLPU[t]);
		printf("rejection probabilities for U, |e_V| = %d, range of values: ", t);
		show_farray(rU[t]);
		printf("%g\n", (n2 - kU + d - t + z) / 3);
	}
}

int main(int argc, char ** argv) {
	struct gengetopt_args_info args_info;

  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1) ;

	init(args_info.params_id_arg);
	show_precomp();
	cleanup();
	cmdline_parser_free(&args_info);
	return 0;
}
