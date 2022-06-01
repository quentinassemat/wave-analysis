#include "wave.h"

/* Wave key i/o */

void wave_keypair_write(wave_sk_t sk, wave_pk_t pk, char * params_id, unsigned long seed);

wave_pk_t wave_publickey(char * params_id, unsigned long key_seed);
wave_sk_t wave_secretkey(char * params_id, unsigned long key_seed);

char * filename_prefix(char * params_id, unsigned long kseed, unsigned long mseed, unsigned long count, char * options);
FILE * fopen_sign(char * prefix);
FILE * fopen_log(char * prefix);
