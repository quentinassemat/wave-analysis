#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cmdline.h"
#include "wave.h"
#include "wave_io.h"

int main(int argc, char ** argv) {
	unsigned long i, j, key_seed, message_seed, fail = 0;
	f3_t * e, * s;
	wave_pk_t pk;
	char * params_id;
	FILE * stream;
	prng_t PRNG;
	struct gengetopt_args_info args_info;

  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1) ;

	verbose = args_info.verbose_arg;
	params_id = args_info.params_id_arg;
	key_seed = args_info.key_arg;

	init_params(params_id);
	pk = wave_publickey(params_id, key_seed);

	if (! args_info.filename_given) {
		printf("verification requires an input file containing signatures\n");
		exit(0);
	}
	stream = fopen(args_info.filename_arg, "r");
	if (stream == NULL) {
		printf("cannot open signature file %s\n", args_info.filename_arg);
		exit(0);
	}

	e = malloc(n);
	for (i = 0; ; ++i) {
		j = fread(&message_seed, sizeof (unsigned long), 1, stream);
		if (j == 0)
			break;
		j = f3_array_read(e, n, stream);
		if (j < n)
			break;
		PRNG = prng_init("message", params_id, message_seed);
		s = f3_randarray(n, PRNG);
		if (! verif(s, e, pk)) {
			printf("verification failed!\n");
			++fail;
		}
		if ((verbose > 2) && i && (i % 10000 == 0))
			printf("%lu\n", i);
		free(s);
		prng_clear(PRNG);
	}
	if (verbose > 0) {
		printf("%lu signature", i);
		if (i > 1)
				printf("s");
		printf(" tested, ");
		if (fail > 0) {
			printf("%lu failure", fail);
			if (fail > 1)
				printf("s");
			printf("\n");
		}
		else
			printf("no failure\n");

	}
	free(e);
	fclose(stream);
	wave_pk_clear(pk);
	cleanup();
	cmdline_parser_free(&args_info);
	return 0;
}
