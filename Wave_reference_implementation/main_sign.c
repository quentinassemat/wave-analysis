#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cmdline.h"
#include "wave.h"
#include "wave_io.h"

void write_log(FILE * stream, unsigned long count) {
	rewind(stream);
	fprintf(stream, "%lu signatures\n", count);
	fprintf(stream, "%lu rejections for V\n%lu rejections for U\n", rejV, rejU);
	fprintf(stream, "%lu Gaussian elimination for V\n%lu Gaussian elimination for U\n", gaussV, gaussU);
	fflush(stream);
}

int main(int argc, char ** argv) {
	unsigned long j, fail = 0, key_seed, message_seed, count;
	f3_t * e, * s;
	wave_sk_t sk;
	wave_pk_t pk;
	char * params_id, * prefix;
	FILE * stream_sign, * stream_log;
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

	sk = wave_secretkey(params_id, key_seed);
	if (args_info.check_flag)
		pk = wave_publickey(params_id, key_seed);

	/* statistics for log file */
	rejV = rejU = gaussU = gaussV = 0;

	prefix = filename_prefix(params_id, key_seed, message_seed, count, NULL);
	stream_sign = (args_info.no_save_flag) ? NULL : fopen_sign(prefix);
	stream_log = (args_info.log_flag) ? fopen_log(prefix) : NULL;
	if (verbose > 0) {
		if (stream_sign != NULL)
			printf("signatures saved in file %s.dat\n", prefix);
		if (stream_log != NULL)
			printf("log file is %s.log\n", prefix);
	}
	for (j = 0; j < count; ++j, ++message_seed) {
		PRNG = prng_init("message", params_id, message_seed);
		s = f3_randarray(n, PRNG);
		prng_clear(PRNG);
		PRNG = prng_init("decode", params_id, message_seed);
		e = sign_wave(s, sk, PRNG);
		prng_clear(PRNG);
		if (args_info.check_flag && (! verif(s, e, pk))) {
			if (verbose > 0)
				printf("verification failed for message seed %lu!\n", message_seed);
			fail++;
		}
		if (stream_sign != NULL) {
			// fwrite(&message_seed, sizeof (unsigned long), 1, stream_sign);
			f3_array_write(e, n, stream_sign);
			fflush(stream_sign);
		}
		if (stream_log != NULL)
			write_log(stream_log, j + 1);
		free(s);
		free(e);
	}
	if (stream_log != NULL) {
		write_log(stream_log, count);
	}
	if (verbose > 0) {
		if (args_info.check_flag) {
			printf("%lu signature", count);
			if (count > 1)
				printf("s");
			printf(" generated, ");
			if (fail > 0) {
				printf("%lu verification", fail);
				if (fail > 1)
					printf("s");
				printf(" failed\n");
			}
			else
				printf("no failure\n");
		}
		if (args_info.stat_given) {
			if (!args_info.check_flag)
				printf("%lu signatures\n", j);
			printf("%lu rejections for V\n%lu rejections for U\n", rejV, rejU);
			printf("%lu Gaussian elimination for V\n%lu Gaussian elimination for U\n", gaussV, gaussU);
		}
	}
	if (stream_sign != NULL)
		fclose(stream_sign);
	if (stream_log != NULL)
		fclose(stream_log);
	free(prefix);
	wave_sk_clear(sk);
	if (args_info.check_flag)
		wave_pk_clear(pk);
	cleanup();
	cmdline_parser_free(&args_info);
	return 0;
}
