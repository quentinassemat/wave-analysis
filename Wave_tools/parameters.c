/**
 * \file parameters.c : files from Wave implementation to use Wave's parameters
 */


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "wave.h"
#include "params.h"

char *current_params_id = NULL;

void init_params(char *params_id)
{
	if ((current_params_id == NULL) || (strcmp(current_params_id, params_id) != 0))
	{
		current_params_id = realloc(current_params_id, strlen(params_id) + 1);
		strcpy(current_params_id, params_id);
	}
	int i = 0;
	while (1)
	{
		if (strcmp(defined_params[i].identifier, "") == 0)
		{
			printf("unknown parameters identifier '%s'\n", params_id);
			exit(0);
		}
		if (strcmp(defined_params[i].identifier, params_id) == 0)
		{
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

void cleanup_params()
{
	if (current_params_id != NULL)
		free(current_params_id);
	current_params_id = NULL;
}
