#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cmdline.h"
#include "wave.h"
#include "wave_io.h"

// en plus
#include "math.h"
#include "time.h"

// #define PRINT_COVARIANCE_MATRIX

#define NB_ANALYSIS_COV 3

int main(int argc, char **argv)
{
	// Prise en compte des arguments/options
	unsigned long i, j, key_seed, verbose, full_test = 0;
	f3_t *e;
	wave_sk_t sk;
	char *params_id;
	FILE *stream;
	struct gengetopt_args_info args_info;

	if (cmdline_parser(argc, argv, &args_info) != 0)
		exit(1);

	params_id = args_info.params_id_arg;
	key_seed = args_info.key_arg;
	verbose = args_info.verbose_arg;
	full_test = args_info.full_test_arg;

	init_params(params_id);

	if (args_info.print_params_flag)
	{
		printf("n, w, kU, kV, d = %d, %d, %d, %d, %d\n", n, w, kU, kV, d);
		exit(0);
	}

	// Start of the analysis

	// Target distribution
	double p0 = ((double)((n - w) * (n - w - 1))) / ((double)(n * (n - 1)));
	double p1 = ((double)(2 * w * (n - w))) / ((double)(n * (n - 1)));
	double p2 = ((double)(w * (w - 1))) / ((double)(n * (n - 1)));

	double proba_class[3] = {p0, p1, p2};
	double proba[9] = {p0, p1 / 4, p1 / 4, p1 / 4, p2 / 4, p2 / 4, p1 / 4, p2 / 4, p2 / 4};
	printf("Target distribution : %f %f %f %f %f %f %f %f %f\n", proba[0], proba[1], proba[2], proba[3], proba[4], proba[5], proba[6], proba[7], proba[8]);

	if (!full_test)
	{
		// order1 is the matrix storing empiric distribution of analysed pairs.
		double **ordre1 = malloc(n * sizeof(double *));

		for (int i1 = 0; i1 < n; i1++)
		{
			ordre1[i1] = malloc(9 * sizeof(double));
			for (int i2 = 0; i2 < 9; i2++)
			{
				ordre1[i1][i2] = 0.0;
			}
		}

		sk = wave_secretkey(params_id, key_seed);

		// Inversion of the secret permutation to detect matched pairs.
		int perm_inv[n];
		for (i = 0; i < n; i++)
		{
			perm_inv[sk.perm[i]] = i;
		}

		if (!args_info.filename_given)
		{
			printf("verification requires an input file containing signatures\n");
			exit(0);
		}
		stream = fopen(args_info.filename_arg, "r");
		if (stream == NULL)
		{
			printf("cannot open signature file %s\n", args_info.filename_arg);
			exit(0);
		}

		// Order 1 analysis

		// Choice of the pairs to analyse i.e. [matched] || [random/not matched]
		time_t seed = time(NULL);
		srand(seed);
		int index_i[n];
		int index_j[n];
		for (int i1 = 0; i1 < n2; i1++)
		{
			// pairs matchées
			index_i[i1] = perm_inv[i1];
			index_j[i1] = perm_inv[i1 + n2];
			// pairs aléatoires (a priori pas matchées)
			index_i[i1 + n2] = rand() % n;
			index_j[i1 + n2] = rand() % n;
			// paires pas aléatoires (pas matchées de manière certaine)
			// index_i[i1 + n2] = perm_inv[i1];
			// index_j[i1 + n2] = perm_inv[(i1 + n2 + 1 >= n) ? (i1 + 1) : (i1 + n2 + 1)];
		}

		printf("Order 1 analysis\n");
		e = malloc(n);
		f3_t ei, ej;
		for (i = 0;; ++i)
		{
			j = f3_array_read(e, n, stream);
			if (j < n)
				break;
			for (int i1 = 0; i1 < n; i1++)
			{
				ei = e[index_i[i1]];
				ej = e[index_j[i1]];
				ordre1[i1][pair_sort(ei, ej)] += 1.0;
			}
		}
		for (int i1 = 0; i1 < n; i1++)
		{
			for (int i2 = 0; i2 < 9; i2++)
			{
				ordre1[i1][i2] /= (double)i;
			}
		}
		fclose(stream);

		// Order 2 analysis

		printf("Order 2 analysis\n");
		// As order 2 analysis implies lots of computation we analyse n2 matched pairs and n2 non-matched pairs.

		stream = fopen(args_info.filename_arg, "r");
		if (stream == NULL)
		{
			printf("cannot open signature file %s\n", args_info.filename_arg);
			exit(0);
		}

		double ****ordre2 = malloc(NB_ANALYSIS_COV * 2 * sizeof(double ***));
		for (int i1 = 0; i1 < NB_ANALYSIS_COV * 2; i1++)
		{
			ordre2[i1] = malloc(NB_ANALYSIS_COV * 2 * sizeof(double **));
			for (int i2 = 0; i2 < NB_ANALYSIS_COV * 2; i2++)
			{
				ordre2[i1][i2] = malloc(9 * sizeof(double *));
				for (int i3 = 0; i3 < 9; i3++)
				{
					ordre2[i1][i2][i3] = malloc(9 * sizeof(double));
					for (int i4 = 0; i4 < 9; i4++)
					{
						ordre2[i1][i2][i3][i4] = 0;
					}
				}
			}
		}

		for (i = 0;; ++i)
		{
			// j = fread(&message_seed, sizeof(unsigned long), 1, stream);
			// if (j == 0)
			// 	break;
			j = f3_array_read(e, n, stream);
			if (j < n)
				break;
			for (int i1 = 0; i1 < NB_ANALYSIS_COV; i1++)
			{
				for (int i2 = 0; i2 < NB_ANALYSIS_COV; i2++)
				{
					ordre2[i1][i2][pair_sort(e[index_i[i1]], e[index_j[i1]])][pair_sort(e[index_i[i2]], e[index_j[i2]])]++;
					ordre2[i1 + NB_ANALYSIS_COV][i2][pair_sort(e[index_i[i1 + n2]], e[index_j[i1 + n2]])][pair_sort(e[index_i[i2]], e[index_j[i2]])]++;
					ordre2[i1][i2 + NB_ANALYSIS_COV][pair_sort(e[index_i[i1]], e[index_j[i1]])][pair_sort(e[index_i[i2 + n2]], e[index_j[i2 + n2]])]++;
					ordre2[i1 + NB_ANALYSIS_COV][i2 + NB_ANALYSIS_COV][pair_sort(e[index_i[i1 + n2]], e[index_j[i1 + n2]])][pair_sort(e[index_i[i2 + n2]], e[index_j[i2 + n2]])]++;
				}
			}
		}

		// ordre2[i][j][k][l] is the number of times the pairs (i,j) took the value (k,l)

		fclose(stream);

		double **covariance = malloc(NB_ANALYSIS_COV * 2 * sizeof(double *));
		for (int i1 = 0; i1 < NB_ANALYSIS_COV * 2; i1++)
		{
			covariance[i1] = malloc(NB_ANALYSIS_COV * 2 * sizeof(double));
			for (int i2 = 0; i2 < NB_ANALYSIS_COV * 2; i2++)
			{
				covariance[i1][i2] = ordre2[i1][i2][0][0];
			}
		}
		// We compute covariance with (e_i, e_j) pairs as Bernoulli variable 1_{(e_i, e_j) = (0, 0)}
		for (int i1 = 0; i1 < NB_ANALYSIS_COV * 2; i1++)
		{
			for (int i2 = 0; i2 < NB_ANALYSIS_COV * 2; i2++)
			{
				covariance[i1][i2] /= (double)i;
			}
		}

		for (int i1 = 0; i1 < NB_ANALYSIS_COV; i1++)
		{
			for (int i2 = 0; i2 < NB_ANALYSIS_COV; i2++)
			{
				covariance[i1][i2] -= ordre1[i1][0] * ordre1[i2][0];
				covariance[i1 + NB_ANALYSIS_COV][i2] -= ordre1[i1 + n2][0] * ordre1[i2][0];
				covariance[i1][i2 + NB_ANALYSIS_COV] -= ordre1[i1][0] * ordre1[i2 + n2][0];
				covariance[i1 + NB_ANALYSIS_COV][i2 + NB_ANALYSIS_COV] -= ordre1[i1 + n2][0] * ordre1[i2 + n2][0];
			}
		}

		// Start of statistical test

		// Choice of parameters:
		// Overall test risk is 0.01 and for that purpose we have (see pdf) (alpha, threshold) = (0.05, 246) (0.5, 2199) (0.95, 4066)

		double alpha = 0.05;
		double eps = sqrt((log(2) - log(alpha)) / (2 * i));
		// attention test_cible n'est pas de niveau alpha

		printf("Bernoulli/Chernoff test:\n");
		printf("alpha : %f, eps : %f\n", alpha, eps);

		int test = test_cible(ordre1, n, proba, eps, 246, verbose);
		printf("test 1 ... %s\n", (test) ? "OK" : "KO");

		test = test_cible_class(ordre1, n, proba_class, eps, 246, verbose);
		printf("test 2 ... %s\n", (test) ? "OK" : "KO");

		test = test_cible_class2(ordre1, n, proba_class, eps, 246, verbose);
		printf("test 3 ... %s\n", (test) ? "OK" : "KO");

		alpha = 0.5;
		eps = sqrt((log(2) - log(alpha)) / (2 * i));
		printf("alpha : %f, eps : %f\n", alpha, eps);

		test = test_cible(ordre1, n, proba, eps, 2199, verbose);
		printf("test 4 ... %s\n", (test) ? "OK" : "KO");

		test = test_cible_class(ordre1, n, proba_class, eps, 2199, verbose);
		printf("test 5 ... %s\n", (test) ? "OK" : "KO");

		test = test_cible_class2(ordre1, n, proba_class, eps, 2199, verbose);
		printf("test 6 ... %s\n", (test) ? "OK" : "KO");

		alpha = 0.95;
		eps = sqrt((log(2) - log(alpha)) / (2 * i));
		printf("alpha : %f, eps : %f\n", alpha, eps);

		test = test_cible(ordre1, n, proba, eps, 4066, verbose);
		printf("test 7 ... %s\n", (test) ? "OK" : "KO");

		test = test_cible_class(ordre1, n, proba_class, eps, 4066, verbose);
		printf("test 8 ... %s\n", (test) ? "OK" : "KO");

		test = test_cible_class2(ordre1, n, proba_class, eps, 4066, verbose);
		printf("test 9 ... %s\n", (test) ? "OK" : "KO");

		printf("\nKhi-Deux test (for %lu signatures):\n", i);
		// 1 - alpha quantial for alpha = 0.95 0.5 0.100 0.050 0.025 0.010 0.005
		double quant[7] = {2.73, 7.34, 13.36, 15.51, 17.53, 20.09, 21.96};

		test = test_adequation(ordre1, n, proba, i, quant[3], 246, verbose);
		printf("test 10 ... %s\n", (test) ? "OK" : "KO");
		test = test_adequation(ordre1, n, proba, i, quant[1], 2199, verbose);
		printf("test 11 ... %s\n", (test) ? "OK" : "KO");
		test = test_adequation(ordre1, n, proba, i, quant[0], 4066, verbose);
		printf("test 12 ... %s\n\n", (test) ? "OK" : "KO");

		printf("Statistical distance test:\n");
		test_distance_stats(ordre1, n, proba, i);

		// Order 2 test
		printf("\nIndependance test :\n");
		test_independance(ordre2, ordre1, NB_ANALYSIS_COV, i);

		printf("Tested pairs:\n");
		printf("matched:");
		for (int i1 = 0; i1 < NB_ANALYSIS_COV; i1++)
		{
			printf("\t(%d, %d)", index_i[i1], index_j[i1]);
		}
		printf("\nrandom/not matched:");
		for (int i1 = 0; i1 < NB_ANALYSIS_COV; i1++)
		{
			printf("\t(%d, %d)", index_i[i1 + n2], index_j[i1 + n2]);
		}
		printf("\n");

#ifdef PRINT_COVARIANCE_MATRIX
		print_covariance(covariance, 2 * NB_ANALYSIS_COV);
#endif

		return 0;
	}
	else
	{
		long unsigned int nb_tot = n*(n - 1) / 2;
		double **order1_tot = malloc(nb_tot * sizeof(double *));
		double **order1_matched = malloc(n2 * sizeof(double *));

		for (int i1 = 0; i1 < nb_tot; i1++)
		{
			order1_tot[i1] = malloc(9 * sizeof(double));
			for (int i2 = 0; i2 < 9; i2++)
			{
				order1_tot[i1][i2] = 0.0;
			}
		}

		for (int i1 = 0; i1 < n2; i1++)
		{
			order1_matched[i1] = malloc(9 * sizeof(double));
			for (int i2 = 0; i2 < 9; i2++)
			{
				order1_matched[i1][i2] = 0.0;
			}
		}

		sk = wave_secretkey(params_id, key_seed);

		// Inversion of the secret permutation to detect matched pairs.
		int perm_inv[n];
		for (i = 0; i < n; i++)
		{
			perm_inv[sk.perm[i]] = i;
		}

		if (!args_info.filename_given)
		{
			printf("verification requires an input file containing signatures\n");
			exit(0);
		}
		stream = fopen(args_info.filename_arg, "r");
		if (stream == NULL)
		{
			printf("cannot open signature file %s\n", args_info.filename_arg);
			exit(0);
		}

		// Order 1 analysis

		// Choice of the pairs to analyse i.e. [matched] || [random/not matched]
		time_t seed = time(NULL);
		srand(seed);
		int index_i[n2];
		int index_j[n2];
		for (int i1 = 0; i1 < n2; i1++)
		{
			// pairs matchées
			index_i[i1] = perm_inv[i1];
			index_j[i1] = perm_inv[i1 + n2];
		}

		printf("Order 1 analysis\n");
		e = malloc(n);
		f3_t ei, ej;
		for (i = 0;; ++i)
		{
			j = f3_array_read(e, n, stream);
			if (j < n)
				break;

			long unsigned int index_tot = 0;

			for (int i1 = 0; i1 < n; i1++)
			{
				for (int i2 = i1 + 1; i2 < n; i2++)
				{
					ei = e[i1];
					ej = e[i2];
					order1_tot[index_tot++][pair_sort(ei, ej)] += 1;
				}
			}
			for (int i1 = 0; i1 < n2; i1 ++) {
				ei = e[index_i[i1]];
				ej = e[index_j[i1]];
				order1_matched[i1][pair_sort(ei, ej)] += 1;
			}
		}
		for (int i1 = 0; i1 < n2; i1++)
		{
			for (int i2 = 0; i2 < 9; i2++)
			{
				order1_matched[i1][i2] /= (double)i;
			}
		}
		for (int i1 = 0; i1 < nb_tot; i1++)
		{
			for (int i2 = 0; i2 < 9; i2++)
			{
				order1_tot[i1][i2] /= (double)i;
			}
		}
		fclose(stream);

		double *t_tot = malloc(nb_tot * sizeof(double));
		double *t_matched = malloc(n2 * sizeof(double));
		for (int i1 = 0; i1 < nb_tot; i1++)
		{
			t_tot[i1] = 0;
		}
		for (int i1 = 0; i1 < n2; i1++)
		{
			t_matched[i1] = 0;
		}
		for (int i1 = 0; i1 < nb_tot; i1++)
		{
			for (int i2 = 0; i2 < 9; i2++)
			{
				t_tot[i1] += pow((order1_tot[i1][i2] - proba[i2]), 2) / proba[i2];
			}
			t_tot[i1] *= i;
		}
		for (int i1 = 0; i1 < n2; i1++)
		{
			for (int i2 = 0; i2 < 9; i2++)
			{
				t_matched[i1] += pow((order1_matched[i1][i2] - proba[i2]), 2) / proba[i2];
			}
			t_matched[i1] *= i;
		}

		FILE* output_matched = NULL;
		FILE* output_tot = NULL;
		output_matched = fopen("../Data/chi2-visualisation/output_matched.txt", "w+");
		output_tot = fopen("../Data/chi2-visualisation/output_tot.txt", "w+");

		for (int i1 = 0; i1 < nb_tot; i1 ++) {
			fprintf(output_tot, "%lf\n", t_tot[i1]);
		}

		for (int i1 = 0; i1 < n2; i1++) {
			fprintf(output_matched, "%lf\n", t_matched[i1]);
		}
		fclose(output_matched);
		fclose(output_tot);

		free(t_tot);
		free(t_matched);
		for (int i1 = 0; i1 < n2; i1 ++) {
			free(order1_matched[i1]);
		}
		for (int i1 = 0; i1 < nb_tot; i1 ++) {
			free(order1_tot[i1]);
		}
		free(order1_tot);
		free(order1_matched);

		return 1;
	}
}