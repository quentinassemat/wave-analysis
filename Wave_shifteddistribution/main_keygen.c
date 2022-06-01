#include <stdlib.h>
#include <stdio.h>
#include "cmdline.h"
#include "wave.h"
#include "wave_io.h"

int main(int argc, char ** argv) {
	unsigned long key_seed;
	char * params_id;
	wave_sk_t sk;
	wave_pk_t pk;
	prng_t PRNG;
	struct gengetopt_args_info args_info;

  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(1) ;

	params_id = args_info.params_id_arg;
	key_seed = args_info.key_arg;
	verbose = args_info.verbose_arg;

	// initialize random generator from seed
	PRNG = prng_init("key", params_id, key_seed);

	init_params(params_id);
	keygen(&sk, &pk, PRNG);
	wave_keypair_write(sk, pk, params_id, key_seed);

	wave_sk_clear(sk);
	wave_pk_clear(pk);
	prng_clear(PRNG);
	cleanup();
	cmdline_parser_free(&args_info);
	return 0;
}
