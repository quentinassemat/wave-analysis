#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cmdline.h"
#include "wave.h"
#include "wave_io.h"


int main(int argc, char ** argv) {
	unsigned long j, key_seed, message_seed, count;
	f3_t * e;
	char * params_id, * prefix;
	FILE * stream_word;
	prng_t PRNG;
	struct gengetopt_args_info args_info;

	if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1);

	verbose = args_info.verbose_arg;
	params_id = args_info.params_id_arg;
	key_seed = args_info.key_arg;
	message_seed = args_info.message_arg;
	count = args_info.number_arg;

	init(params_id);

	if (args_info.print_params_flag) {
		printf("n, w, kU, kV, d = %d, %d, %d, %d, %d\n", n, w, kU, kV, d);
		exit(0);
	}


	prefix = filename_prefix(params_id, key_seed, message_seed, count, NULL);
	stream_word = (args_info.no_save_flag) ? NULL : fopen_sign(prefix);

	if (verbose > 0) {
		if (stream_word != NULL)
			printf("signatures saved in file %s.dat\n", prefix);
	}
	for (j = 0; j < count; ++j, ++message_seed) {
		PRNG = prng_init("mot", params_id, message_seed);
		e = rand_word(n, w, PRNG);
		prng_clear(PRNG);
		if (stream_word != NULL) {
			// fwrite(&message_seed, sizeof (unsigned long), 1, stream_word);
			f3_array_write(e, n, stream_word);
			fflush(stream_word);
		}
		free(e);
	}
	if (stream_word != NULL)
		fclose(stream_word);
	free(prefix);
	cleanup();
	cmdline_parser_free(&args_info);
	return 0;
}
