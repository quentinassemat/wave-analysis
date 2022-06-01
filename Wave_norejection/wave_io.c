#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "wave.h"

char * pk_filename(char * params_id, unsigned long seed) {
	char * filename;

#define PK_FILE_FORMAT "../Data/%s/wave_pk%s_%lu.dat"
	// 20 is the decimal size of the largest 'unsigned long'
	filename = malloc(strlen(PK_FILE_FORMAT) + 2 * strlen(params_id) + 20);
	sprintf(filename, PK_FILE_FORMAT, params_id, params_id, seed);
	return filename;
}

wave_pk_t wave_publickey_read(char * filename) {
	FILE * stream = fopen(filename, "r");
	if (stream == NULL) {
		printf("Cannot open file %s!\n", filename);
		exit(0);
	}
	int j = 0;
	wave_pk_t pk;

	j += fread_publickey(&pk, stream);
	fclose(stream);

	return pk;
}

int wave_publickey_write(wave_pk_t pk, char * filename) {
	FILE * stream = fopen(filename, "w");
	if (stream == NULL) {
		printf("Cannot open file %s!\n", filename);
		exit(0);
	}
	int j = 0;
	j += fwrite_publickey(&pk, stream);
	fclose(stream);
	if (verbose > 0)
		printf("public key saved in file %s (%d bytes)\n", filename, j);
	return j;
}

wave_pk_t wave_publickey(char * params_id, unsigned long key_seed) {
	char * filename = pk_filename(params_id, key_seed);
	wave_pk_t pk = wave_publickey_read(filename);
	free(filename);
	return pk;
}

char * sk_filename(char * params_id, unsigned long seed) {
	char * filename;

#define SK_FILE_FORMAT "../Data/%s/wave_sk%s_%lu.dat"
	filename = malloc(strlen(SK_FILE_FORMAT) + 2 * strlen(params_id) + 11);
	sprintf(filename, SK_FILE_FORMAT, params_id, params_id, seed);
	return filename;
}

int wave_secretkey_write(wave_sk_t sk, char * filename) {
	FILE * stream = fopen(filename, "w");
	if (stream == NULL) {
		printf("Cannot open file %s!\n", filename);
		exit(0);
	}
	int j = 0;
	j += fwrite_secretkey(&sk, stream);
	fclose(stream);
	if (verbose > 0)
		printf("secret key saved in file %s (%d bytes)\n", filename, j);
	return j;
}

wave_sk_t wave_secretkey_read(char * filename) {
	FILE * stream = fopen(filename, "r");
	if (stream == NULL) {
		printf("Cannot open file %s!\n", filename);
		exit(0);
	}
	int j = 0;
	wave_sk_t sk;
	j += fread_secretkey(&sk, stream);
	fclose(stream);
	return sk;
}

wave_sk_t wave_secretkey(char * params_id, unsigned long key_seed) {
	char * filename = sk_filename(params_id, key_seed);
	wave_sk_t sk = wave_secretkey_read(filename);
	free(filename);
	return sk;
}

void wave_keypair_write(wave_sk_t sk, wave_pk_t pk, char * params_id, unsigned long seed) {
	char * filename;

	filename = sk_filename(params_id, seed);
	wave_secretkey_write(sk, filename);
	free(filename);

	filename = pk_filename(params_id, seed);
	wave_publickey_write(pk, filename);
	free(filename);
}

char * filename_prefix(char * params_id, unsigned long kseed, unsigned long mseed, unsigned long count, char * options) {
	char * filename;
#define FILENAME_FORMAT "../Data/%s/sign_norejection%s_%s_%lu_%lu-%lu"
	// 20 is the decimal size of the largest 'unsigned long'
	if ((options == NULL) || (strcmp(options, "0") == 0))
		options = "";
	filename = malloc(strlen(FILENAME_FORMAT) + strlen(options) + 2 * strlen(params_id) + 3 * 20);
	sprintf(filename, FILENAME_FORMAT, params_id, options, params_id, kseed, mseed, mseed + count - 1);
	return filename;
}

FILE * fopen_aux(char * prefix, char * suffix) {
	FILE * stream;
	char * filename;
	filename = malloc(strlen(prefix) + strlen(suffix) + 1);
	sprintf(filename, "%s%s", prefix, suffix);
	stream = fopen(filename, "w");
	if (stream == NULL) {
		printf("cannot open file %s\n", filename);
		exit(0);
	}
	free(filename);

	return stream;
}

FILE * fopen_sign(char * prefix) {
	return fopen_aux(prefix, ".dat");
}

FILE * fopen_log(char * prefix) {
	return fopen_aux(prefix, ".log");
}
