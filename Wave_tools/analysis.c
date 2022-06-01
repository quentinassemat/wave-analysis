/**
 * \file analysis.c : files which contains all the different statistical tests
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cmdline.h"
#include "wave.h"
#include "wave_io.h"
#include "math.h"

// Order 1 tests :

/**
 * @brief test if the distribution (notably among the matched pairs) is the expected one.
 * @param order1 : matrix storing empiric distribution of analysed pairs.
 * @param lenght the dimension of order1
 * @param proba the target distribution
 * @param eps the epsilon choose to set the overall test to the target risk.
 * @param threshold the threshold to set the overall test to the target risk.
 * @param verbose to display more details
 * @return return 0 if the test fail and 1 otherwise
 */
int test_cible(double **ordre1, int lenght, double *proba, double eps, int threshold, int verbose)
{
	int anomalies_matched[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // anomalies[i] est le nombre de rejet pour le test en étudiant la variable de bernoulli 1_{x, y = a, b}
	int anomalies_rand[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};	// anomalies[i] est le nombre de rejet pour le test en étudiant la variable de bernoulli 1_{x, y = a, b}
	for (int i1 = 0; i1 < 9; i1++)
	{
		anomalies_matched[i1] = 0;
		anomalies_rand[i1] = 0;
	}
	for (int i1 = 0; i1 < lenght / 2; i1++)
	{
		for (int i2 = 0; i2 < 9; i2++)
		{
			int suspect = ((ordre1[i1][i2] > proba[i2] + eps) || (ordre1[i1][i2] < proba[i2] - eps));
			anomalies_matched[i2] += suspect;
			suspect = ((ordre1[i1 + lenght / 2][i2] > proba[i2] + eps) || (ordre1[i1 + lenght / 2][i2] < proba[i2] - eps));
			anomalies_rand[i2] += suspect;
		}
	}
	if (verbose > 0)
	{
		printf("Rejection among matched: %d %d %d %d %d %d %d %d %d\n", anomalies_matched[0], anomalies_matched[1], anomalies_matched[2], anomalies_matched[3], anomalies_matched[4], anomalies_matched[5], anomalies_matched[6], anomalies_matched[7], anomalies_matched[8]);
		printf("Rejection among random/non matched: %d %d %d %d %d %d %d %d %d\n", anomalies_rand[0], anomalies_rand[1], anomalies_rand[2], anomalies_rand[3], anomalies_rand[4], anomalies_rand[5], anomalies_rand[6], anomalies_rand[7], anomalies_rand[8]);
	}
	int fail = 1;
	for (int i1 = 0; i1 < 9; i1++)
	{
		if (anomalies_matched[i1] > threshold)
		{
			fail = 0;
		}
		if (anomalies_rand[i1] > threshold)
		{
			fail = 0;
		}
	}
	return fail;
}

/**
 * @brief test if the distribution (notably among the matched pairs) is the expected one. This test gathered possible outcomes of pairs by class (see pdf).
 * @param order1 : matrix storing empiric distribution of analysed pairs.
 * @param lenght the dimension of order1
 * @param proba_class the target distribution
 * @param eps the epsilon choosed to set the overall test to the target risk.
 * @param threshold the threshold choosed to set the overall test to the target risk.
 * @param verbose to display more details
 * @return return 0 if the test fail and 1 otherwise
 */
int test_cible_class(double **ordre1, int lenght, double *proba_class, double eps, int threshold, int verbose)
{
	int anomalies_matched[3] = {0, 0, 0}; // anomalies[i] est le nombre de rejet pour le test en étudiant la variable de bernoulli 1_{x, y = a, b}
	int anomalies_rand[3] = {0, 0, 0}; // anomalies[i] est le nombre de rejet pour le test en étudiant la variable de bernoulli 1_{x, y = a, b}
	for (int i1 = 0; i1 < 3; i1++)
	{
		anomalies_matched[i1] = 0;
		anomalies_rand[i1] = 0;
	}
	for (int i1 = 0; i1 < lenght / 2; i1++)
	{
		double obs = ordre1[i1][0];
		int suspect = ((obs > proba_class[0] + eps) || (obs < proba_class[0] - eps));
		anomalies_matched[0] += suspect;
		obs = ordre1[i1][1] + ordre1[i1][2] + ordre1[i1][3] + ordre1[i1][6];
		suspect = ((obs > proba_class[1] + eps) || (obs < proba_class[1] - eps));
		anomalies_matched[1] += suspect;
		obs = ordre1[i1][4] + ordre1[i1][5] + ordre1[i1][7] + ordre1[i1][8];
		suspect = ((obs > proba_class[2] + eps) || (obs < proba_class[2] - eps));
		anomalies_matched[2] += suspect;

		obs = ordre1[i1 + lenght / 2][0];
		suspect = ((obs > proba_class[0] + eps) || (obs < proba_class[0] - eps));
		anomalies_rand[0] += suspect;
		obs = ordre1[i1 + lenght / 2][1] + ordre1[i1 + lenght / 2][2] + ordre1[i1 + lenght / 2][3] + ordre1[i1 + lenght / 2][6];
		suspect = ((obs > proba_class[1] + eps) || (obs < proba_class[1] - eps));
		anomalies_rand[1] += suspect;
		obs = ordre1[i1 + lenght / 2][4] + ordre1[i1 + lenght / 2][5] + ordre1[i1 + lenght / 2][7] + ordre1[i1 + lenght / 2][8];
		suspect = ((obs > proba_class[2] + eps) || (obs < proba_class[2] - eps));
		anomalies_rand[2] += suspect;
	}
	if (verbose > 0)
	{
		printf("Rejection among matched: %d %d %d\n", anomalies_matched[0], anomalies_matched[1], anomalies_matched[2]);
		printf("Rejection among random/non matched: %d %d %d\n", anomalies_rand[0], anomalies_rand[1], anomalies_rand[2]);
	}
	int fail = 1;
	for (int i1 = 0; i1 < 3; i1++)
	{
		if (anomalies_matched[i1] > threshold)
		{
			fail = 0;
		}
		if (anomalies_rand[i1] > threshold)
		{
			fail = 0;
		}
	}
	return fail;
}

/**
 * @brief test if the distribution (notably among the matched pairs) is the expected one. This test gathered possible outcomes of pairs by another class partition (see pdf).
 * @param order1 : matrix storing empiric distribution of analysed pairs.
 * @param lenght the dimension of order1
 * @param proba_class the target distribution
 * @param eps the epsilon choosed to set the overall test to the target risk.
 * @param threshold the threshold choosed to set the overall test to the target risk.
 * @param verbose to display more details
 * @return return 0 if the test fail and 1 otherwise
 */
int test_cible_class2(double **ordre1, int lenght, double *proba_class, double eps, int threshold, int verbose)
{
	int anomalies_matched[2] = {0, 0}; // anomalies[i] est le nombre de rejet pour le test en étudiant la variable de bernoulli 1_{x, y = a, b}
	int anomalies_rand[2] = {0, 0}; // anomalies[i] est le nombre de rejet pour le test en étudiant la variable de bernoulli 1_{x, y = a, b}
	for (int i1 = 0; i1 < 2; i1++)
	{
		anomalies_matched[i1] = 0;
		anomalies_rand[i1] = 0;
	}
	for (int i1 = 0; i1 < lenght / 2; i1++)
	{
		double obs = ordre1[i1][5] + ordre1[i1][7];
		int suspect = ((obs > 0.5 * proba_class[2] + eps) || (obs < 0.5 * proba_class[2] - eps));
		anomalies_matched[0] += suspect;
		obs = ordre1[i1][4] + ordre1[i1][8];
		suspect = ((obs > 0.5 * proba_class[2] + eps) || (obs < 0.5 * proba_class[2] - eps));
		anomalies_matched[1] += suspect;

		obs = ordre1[i1 + lenght / 2][5] + ordre1[i1 + lenght / 2][7];
		suspect = ((obs > 0.5 * proba_class[2] + eps) || (obs < 0.5 * proba_class[2] - eps));
		anomalies_rand[0] += suspect;
		obs = ordre1[i1 + lenght / 2][4] + ordre1[i1 + lenght / 2][8];
		suspect = ((obs > 0.5 * proba_class[2] + eps) || (obs < 0.5 * proba_class[2] - eps));
		anomalies_rand[1] += suspect;
	}
	if (verbose > 0)
	{
		printf("Rejection among matched: %d %d\n", anomalies_matched[0], anomalies_matched[1]);
		printf("Rejection among random/non matched: %d %d\n", anomalies_rand[0], anomalies_rand[1]);
	}
	int fail = 1;
	for (int i1 = 0; i1 < 2; i1++)
	{
		if (anomalies_matched[i1] > threshold)
		{
			fail = 0;
		}
		if (anomalies_rand[i1] > threshold)
		{
			fail = 0;
		}
	}
	return fail;
}

// TEST KHI-SQUARE

/**
 * @brief test if the distribution (notably among the matched pairs) is the expected one by performing a Khi-Square test. 
 * @param order1 : matrix storing empiric distribution of analysed pairs.
 * @param lenght the dimension of order1
 * @param proba the target distribution
 * @param nb_signature the number of signature tested (necessary to compute the khi-square statistic)
 * @param quantile the quantile choosed to set the overall test to the target risk.
 * @param threshold choosed to set the overall test to the target risk.
 * @param verbose to display more details
 * @return return 0 if the test fail and 1 otherwise
 */
int test_adequation(double **ordre1, int lenght, double *proba, int nb_signature, double quantile, int threshold, int verbose)
{
	double *t = malloc(lenght * sizeof(double));
	for (int i1 = 0; i1 < lenght; i1++)
	{
		t[i1] = 0;
	}
	for (int i1 = 0; i1 < lenght; i1++)
	{
		for (int i2 = 0; i2 < 9; i2++)
		{
			t[i1] += pow((ordre1[i1][i2] - proba[i2]), 2) / proba[i2];
		}
		t[i1] *= nb_signature;
	}
	int rejet_matched = 0;
	int rejet_rand = 0;
	for (int i1 = 0; i1 < lenght / 2; i1++)
	{
		if (t[i1] > quantile)
		{
			rejet_matched++;
		}
		if (t[i1 + (lenght / 2)] > quantile)
		{
			rejet_rand++;
		}
	}
	int fail = 1;
	if (rejet_matched > threshold || rejet_rand > threshold)
	{
		fail = 0;
	}
	if (verbose > 0)
	{
		printf("Rejection among matched: %d\n", rejet_matched);
		printf("Rejection among random/non matched: %d\n", rejet_rand);
	}
	return fail;
}

// TEST DISTANCE STATISTIQUES

/**
 * @brief test if the distribution (notably among the matched pairs) is the expected one by studying various statistical distance (see pdf)
 * @param order1 : matrix storing empiric distribution of analysed pairs.
 * @param lenght the dimension of order1
 * @param proba the target distribution
 * @return return 0 if the test fail and 1 otherwise
 */
int test_distance_stats(double **ordre1, int lenght, double *proba, int nb_signature)
{
	// Différentes distances statistiques :

	double *dist_kl = malloc(lenght * sizeof(double));	 // somme des p log (p/p_cible) (distnace par rapport à la cible)
	double *dist_stat = malloc(lenght * sizeof(double)); // sommes valeurs absolues des différences entre les distributions
	double *dist_hell = malloc(lenght * sizeof(double)); // distnace euclidienne des racines des vecteurs de probabilités

	for (int i1 = 0; i1 < lenght; i1++)
	{
		dist_kl[i1] = 0;
		dist_stat[i1] = 0;
		dist_hell[i1] = 0;
	}

	for (int i1 = 0; i1 < lenght; i1++)
	{
		for (int i2 = 0; i2 < 9; i2++)
		{
			dist_stat[i1] += fabs(ordre1[i1][i2] - proba[i2]);
			dist_kl[i1] += ordre1[i1][i2] * log(ordre1[i1][i2] / proba[i2]);
			dist_hell[i1] += pow(sqrt(ordre1[i1][i2]) - sqrt(proba[i2]), 2);
		}
		dist_hell[i1] = (1 / sqrt(2)) * sqrt(dist_hell[i1]);
	}
	double dist_kl_mean_matched = 0, dist_hell_mean_matched = 0, dist_stat_mean_matched = 0;
	double dist_kl_mean_rand = 0, dist_hell_mean_rand = 0, dist_stat_mean_rand = 0;
	double dist_kl_sup_matched = 0, dist_hell_sup_matched = 0, dist_stat_sup_matched = 0;
	double dist_kl_sup_rand = 0, dist_hell_sup_rand = 0, dist_stat_sup_rand = 0;

	for (int i1 = 0; i1 < lenght / 2; i1++)
	{
		dist_kl_mean_matched += dist_kl[i1];
		dist_kl_mean_rand += dist_kl[i1 + (lenght / 2)];
		dist_hell_mean_matched += dist_hell[i1];
		dist_hell_mean_rand += dist_hell[i1 + (lenght / 2)];
		dist_stat_mean_matched += dist_stat[i1];
		dist_stat_mean_rand += dist_stat[i1 + (lenght / 2)];
		if (dist_kl[i1] > dist_kl_sup_matched)
		{
			dist_kl_sup_matched = dist_kl[i1];
		}
		if (dist_kl[i1 + (lenght / 2)] > dist_kl_sup_rand)
		{
			dist_kl_sup_rand = dist_kl[i1 + (lenght / 2)];
		}
		if (dist_hell[i1] > dist_hell_sup_matched)
		{
			dist_hell_sup_matched = dist_hell[i1];
		}
		if (dist_hell[i1 + (lenght / 2)] > dist_hell_sup_rand)
		{
			dist_hell_sup_rand = dist_hell[i1 + (lenght / 2)];
		}
		if (dist_stat[i1] > dist_stat_sup_matched)
		{
			dist_stat_sup_matched = dist_stat[i1];
		}
		if (dist_stat[i1 + (lenght / 2)] > dist_stat_sup_rand)
		{
			dist_stat_sup_rand = dist_stat[i1 + (lenght / 2)];
		}
	}
	dist_kl_mean_matched /= lenght / 2;
	dist_hell_mean_matched /= lenght / 2;
	dist_stat_mean_matched /= lenght / 2;
	dist_kl_mean_rand /= lenght / 2;
	dist_hell_mean_rand /= lenght / 2;
	dist_stat_mean_rand /= lenght / 2;
	printf("Mean distance among matched pairs: kl : %f    hell : %f    stat : %f\n", dist_kl_mean_matched, dist_hell_mean_matched, dist_stat_mean_matched);
	printf("Mean distance among random/non matched pairs: kl : %f    hell : %f    stat : %f\n\n", dist_kl_mean_rand, dist_hell_mean_rand, dist_stat_mean_rand);

	printf("Supp distance among matched pairs  kl : %f    hell : %f    stat : %f\n", dist_kl_sup_matched, dist_hell_sup_matched, dist_stat_sup_matched);
	printf("Supp distance pairs among random/non matched pairs: kl : %f    hell : %f    stat : %f\n", dist_kl_sup_rand, dist_hell_sup_rand, dist_stat_sup_rand);
	return 0;
}

// Order 2 test

/**
 * @brief test if pairs are independant by performing a Khi-Square test. 
 * @param order1 : matrix storing empiric distribution of analysed pairs.
 * @param order2 : matrix order2 analysis
 * @param nb_cov the number of pairs analysed with order 2.
 * @param nb_signature the number of signature tested (necessary to compute the khi-square statistic)
 * @return return 0 if the test fail and 1 otherwise
 */
int test_independance(double ****ordre2, double **ordre1, int nb_cov, int nb_signatures)
{
	double quant[5] = {78.8596, 83.6753, 88.0041, 93.2169, 96.8781}; // quantile 1 - alpha pour alpha = 0.100 0.050 0.025 0.010 0.005
	int nb_val = 0;
	int nb_rejet = 0;
	int nb_rejet_matched = 0;
	int nb_rejet_rand = 0;
	int nb_rejet_matched_rand = 0;
	double **t = malloc((nb_cov * 2) * sizeof(double *));
	for (int i1 = 0; i1 < nb_cov * 2 - 1; i1++)
	{
		t[i1] = malloc((nb_cov * 2 - 1 - i1) * sizeof(double));
		for (int i2 = 0; i2 < nb_cov * 2 - 1 - i1; i2++)
		{
			nb_val += 1;
			t[i1][i2] = 0;
		}
	}
	// en t[i1][i2] est stocké t_{i1, i2 + i1 + 1} car on ne test pas l'indépendance avec soi même
	for (int i1 = 0; i1 < nb_cov * 2 - 1; i1++)
	{
		for (int i2 = 0; i2 < nb_cov * 2 - i1 - 1; i2++)
		{
			for (int i3 = 0; i3 < 9; i3++)
			{
				for (int i4 = 0; i4 < 9; i4++)
				{
					int j1 = (i1 >= nb_cov) ? i1 + n2 - nb_cov : i1;							//  (i1 >= nb_cov) ? i1 + n2 - nb_cov : i1;
					int j2 = (i2 + i1 + 1 >= nb_cov) ? i2 + i1 + 1 + n2 - nb_cov : i2 + i1 + 1; // (i2 >= nb_cov) ? i2 + i1 + 1 + n2 - nb_cov : i2 + i1 + 1;
					double eij = ordre1[j1][i3] * ordre1[j2][i4] * nb_signatures;
					t[i1][i2] += pow(ordre2[i1][i2 + i1 + 1][i3][i4] - eij, 2) / eij;
				}
			}
		}
	}
	for (int i1 = 0; i1 < nb_cov * 2 - 1; i1++)
	{
#ifdef PRINT_KHIDEUX_IND
		for (int i2 = 0; i2 < i1; i2++)
		{
			printf("       ");
		}
#endif
		for (int i2 = 0; i2 < nb_cov * 2 - i1 - 1; i2++)
		{
			if (t[i1][i2] > quant[4])
			{
				nb_rejet += 1;
				if (i1 + i2 + 1 < nb_cov)
				{
					nb_rejet_matched += 1;
				}
				else
				{
					if (i1 >= nb_cov)
					{
						nb_rejet_rand += 1;
					}
					else
					{
						nb_rejet_matched_rand += 1;
					}
				}
			}
#ifdef PRINT_KHIDEUX_IND
			(t[i1][i2] >= 100) ? printf(" %.2f", t[i1][i2]) : printf("  %.2f", t[i1][i2]);
#endif
		}
#ifdef PRINT_KHIDEUX_IND
		printf("\n");
#endif
	}
	int nb_tot_matched = (nb_cov * (nb_cov - 1)) / 2;
	printf("Rejection : %f\n", ((double)nb_rejet / (double)nb_val) * 100);
	printf("Rejection rejet_matched : %f\n", ((double)nb_rejet_matched / (double)nb_tot_matched) * 100);
	printf("Rejection rejet_rand : %f\n", ((double)nb_rejet_rand / (double)nb_tot_matched) * 100);
	printf("Rejection rejet_matched_rand : %f\n", ((double)nb_rejet_matched_rand / (double)(nb_cov * nb_cov)) * 100);
	return 0;
}

void print_covariance(double **covariance, int lenght)
{
	for (int i = 0; i < lenght; i++)
	{
		for (int j = 0; j < lenght; j++)
		{
			printf("%f ", covariance[i][j]);
		}
		printf("\n");
	}
}

// HELPERS
int pair_sort(f3_t x, f3_t y)
{
	if f3_iszero (x)
	{
		if f3_iszero (y)
		{
			return 0;
		}
		if f3_isone (y)
		{
			return 1;
		}
		else
		{
			return 2;
		}
	}
	if f3_isone (x)
	{
		if f3_iszero (y)
		{
			return 3;
		}
		if f3_isone (y)
		{
			return 4;
		}
		else
		{
			return 5;
		}
	}
	else
	{
		if f3_iszero (y)
		{
			return 6;
		}
		if f3_isone (y)
		{
			return 7;
		}
		else
		{
			return 8;
		}
	}
}

int pair_match(f3_t x, f3_t y)
{
	if f3_iszero (x)
	{
		if f3_iszero (y)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
	if f3_isone (x)
	{
		if f3_iszero (y)
		{
			return 1;
		}
		else
		{
			return 2;
		}
	}
	else
	{
		if f3_iszero (y)
		{
			return 1;
		}
		else
		{
			return 2;
		}
	}
}

// Pas de fond théorique précis mais en même temps utile : à voir

// int test_comp(double **ordre1, int lenght, double *proba, double eps, double threshold, int verbose)
// {
// 	int *anomalies = malloc(lenght * sizeof(int));
// 	for (int i1 = 0; i1 < lenght; i1++)
// 	{
// 		anomalies[i1] = 0;
// 	}
// 	for (int i1 = 0; i1 < lenght; i1++)
// 	{
// 		for (int i2 = 0; i2 < 9; i2++)
// 		{
// 			int suspect = ((ordre1[i1][i2] > proba[i2] + eps) || (ordre1[i1][i2] < proba[i2] - eps));
// 			anomalies[i1] += suspect;
// 		}
// 	}
// 	int stat_anomalies_matched[10];
// 	int stat_anomalies_rand[10];
// 	for (int i1 = 0; i1 < 10; i1++)
// 	{
// 		stat_anomalies_matched[i1] = 0;
// 		stat_anomalies_rand[i1] = 0;
// 	}
// 	for (int i1 = 0; i1 < lenght / 2; i1++)
// 	{
// 		stat_anomalies_matched[anomalies[i1]]++;
// 		stat_anomalies_rand[anomalies[i1 + (lenght / 2)]]++;
// 	}
// 	if (verbose > 0)
// 	{
// 		printf("matched 0 : %d, 1 : %d, 2 : %d, 3 : %d, 4 : %d, 5 : %d, 6 : %d, 7 : %d, 8 : %d, 9 : %d\n", stat_anomalies_matched[0], stat_anomalies_matched[1], stat_anomalies_matched[2], stat_anomalies_matched[3], stat_anomalies_matched[4], stat_anomalies_matched[5], stat_anomalies_matched[6], stat_anomalies_matched[7], stat_anomalies_matched[8], stat_anomalies_matched[9]);
// 		printf("rand    0 : %d, 1 : %d, 2 : %d, 3 : %d, 4 : %d, 5 : %d, 6 : %d, 7 : %d, 8 : %d, 9 : %d\n", stat_anomalies_rand[0], stat_anomalies_rand[1], stat_anomalies_rand[2], stat_anomalies_rand[3], stat_anomalies_rand[4], stat_anomalies_rand[5], stat_anomalies_rand[6], stat_anomalies_rand[7], stat_anomalies_rand[8], stat_anomalies_rand[9]);
// 	}
// 	int nb_anomalies_matched = lenght - stat_anomalies_matched[0];
// 	int nb_anomalies_rand = lenght - stat_anomalies_rand[0];
// 	return (nb_anomalies_matched > threshold * nb_anomalies_rand) ? 0 : 1;
// }

// int test_comp_class(double **ordre1, int lenght, double *proba_class, double eps, double threshold, int verbose)
// {
// 	int *anomalies = malloc(lenght * sizeof(int));
// 	for (int i1 = 0; i1 < lenght; i1++)
// 	{
// 		anomalies[i1] = 0;
// 	}
// 	for (int i1 = 0; i1 < lenght; i1++)
// 	{
// 		double obs = ordre1[i1][0];
// 		int suspect = ((obs > proba_class[0] + eps) || (obs < proba_class[0] - eps));
// 		anomalies[i1] += suspect;
// 		obs = ordre1[i1][1] + ordre1[i1][2] + ordre1[i1][3] + ordre1[i1][6];
// 		suspect = ((obs > proba_class[1] + eps) || (obs < proba_class[1] - eps));
// 		anomalies[i1] += suspect;
// 		obs = ordre1[i1][4] + ordre1[i1][5] + ordre1[i1][7] + ordre1[i1][8];
// 		suspect = ((obs > proba_class[2] + eps) || (obs < proba_class[2] - eps));
// 		anomalies[i1] += suspect;
// 	}
// 	int stat_anomalies_matched[4];
// 	int stat_anomalies_rand[4];
// 	for (int i1 = 0; i1 < 4; i1++)
// 	{
// 		stat_anomalies_matched[i1] = 0;
// 		stat_anomalies_rand[i1] = 0;
// 	}
// 	for (int i1 = 0; i1 < lenght / 2; i1++)
// 	{
// 		stat_anomalies_matched[anomalies[i1]]++;
// 		stat_anomalies_rand[anomalies[i1 + (lenght / 2)]]++;
// 	}
// 	if (verbose > 0)
// 	{
// 		printf("matched 0 : %d, 1 : %d, 2 : %d, 3 : %d\n", stat_anomalies_matched[0], stat_anomalies_matched[1], stat_anomalies_matched[2], stat_anomalies_matched[3]);
// 		printf("rand    0 : %d, 1 : %d, 2 : %d, 3 : %d\n", stat_anomalies_rand[0], stat_anomalies_rand[1], stat_anomalies_rand[2], stat_anomalies_rand[3]);
// 	}
// 	int nb_anomalies_matched = lenght - stat_anomalies_matched[0];
// 	int nb_anomalies_rand = lenght - stat_anomalies_rand[0];
// 	return (nb_anomalies_matched > threshold * nb_anomalies_rand) ? 0 : 1;
// }

// int test_comp_class2(double **ordre1, int lenght, double *proba_class, double eps, double threshold, int verbose)
// {
// 	int *anomalies = malloc(lenght * sizeof(int));
// 	for (int i1 = 0; i1 < lenght; i1++)
// 	{
// 		anomalies[i1] = 0;
// 	}
// 	for (int i1 = 0; i1 < lenght; i1++)
// 	{
// 		double obs = ordre1[i1][5] + ordre1[i1][7]; // (1, 2) ou (2, 1)
// 		int suspect = ((obs > 0.5 * proba_class[2] + eps) || (obs < 0.5 * proba_class[2] - eps));
// 		anomalies[i1] += suspect;
// 		obs = ordre1[i1][4] + ordre1[i1][8]; // (1, 1) ou (2, 2)
// 		suspect = ((obs > 0.5 * proba_class[2] + eps) || (obs < 0.5 * proba_class[2] - eps));
// 		anomalies[i1] += suspect;
// 	}
// 	int stat_anomalies_matched[3];
// 	int stat_anomalies_rand[3];
// 	for (int i1 = 0; i1 < 3; i1++)
// 	{
// 		stat_anomalies_matched[i1] = 0;
// 		stat_anomalies_rand[i1] = 0;
// 	}
// 	for (int i1 = 0; i1 < lenght / 2; i1++)
// 	{
// 		stat_anomalies_matched[anomalies[i1]]++;
// 		stat_anomalies_rand[anomalies[i1 + (lenght / 2)]]++;
// 	}
// 	if (verbose > 0)
// 	{
// 		printf("matched 0 : %d, 1 : %d, 2 : %d\n", stat_anomalies_matched[0], stat_anomalies_matched[1], stat_anomalies_matched[2]);
// 		printf("rand    0 : %d, 1 : %d, 2 : %d\n", stat_anomalies_rand[0], stat_anomalies_rand[1], stat_anomalies_rand[2]);
// 	}
// 	int nb_anomalies_matched = lenght - stat_anomalies_matched[0];
// 	int nb_anomalies_rand = lenght - stat_anomalies_rand[0];
// 	return (nb_anomalies_matched > threshold * nb_anomalies_rand) ? 0 : 1;
// }