#include <stdlib.h>
#include "wave.h"

int verif(f3_t * s, f3_t * e, mf3_t pk) {
	if (f3_array_weight(e, n) != w)
		return 0;
	f3_t * x = mf3_ma_mul(pk, e + (n - k));
	for (int i = 0; i < n - k; ++i) {
		x[i] = f3_add(x[i], e[i]);
		if (x[i] != s[i]) {
			free(x);
			return 0;
		}
	}
	free(x);
	return 1;
}
