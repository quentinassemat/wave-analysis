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

char * filename_prefix(char * params_id, unsigned long kseed, unsigned long mseed, unsigned long count, char * options) {
	char * filename;
#define FILENAME_FORMAT "../Data/%s/word%s_%s_%lu_%lu-%lu"
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
