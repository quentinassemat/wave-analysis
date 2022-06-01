#include "wave.h"

/* identifier, block length, error weight, dimension of U, dimension of V, gap */
wave_params_t defined_params[] = (wave_params_t[]){
	{"128g", 8492, 7980, 3558, 2047, 81},
	{"96g", 6368, 5984, 2668, 1535, 61},
	{"80g", 5308, 4988, 2224, 1280, 50},
	{"64g", 4246, 3990, 1779, 1024, 40},
	{"", 0, 0, 0, 0, 0} // add new parameters sets above and keep this sentinel row last
};
