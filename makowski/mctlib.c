#include <math.h>
#include <omp.h>

#define CHUNKSIZE 1000

// ------------------------------------------//
// ------------------------------------------//

void calculate_vectors(double *var0, double *var1, double *var2, double *var3, double *var4, long int n, double *var5, double *var6, double *var7, double *var8, double *var9)
{
    int chunk = CHUNKSIZE;

    register long int i, j;
    register double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;

#pragma omp for schedule(dynamic, chunk)
    for (i = 0; i < n; i++)
    {
        tmp0 = var0[i];
        tmp1 = var1[i];
        tmp2 = var2[i];
        tmp3 = var3[i];

        tmp5 = sin(tmp1) + cos(tmp0);
        tmp6 = exp(tmp2) + (1 / exp(var4[i]));
        tmp7 = sin(tmp1) * cos(tmp0) + cos(tmp3) * sin(tmp2);
        tmp8 = hypot(tmp2, tmp1);
        tmp9 = cbrt(tmp3);

        var5[i] = tmp5;
        var6[i] = tmp6;
        var7[i] = tmp7;
        var8[i] = tmp8;
        var9[i] = tmp9;
    }
}

// ------------------------------------------//
// ------------------------------------------//

double calculate_mean(double *var, long int n)
{
    int chunk = CHUNKSIZE;

    register long int i;

    double mean_sum = 0.;

#pragma omp for schedule(dynamic, chunk)
    for (i = 0; i < n; i++)
        mean_sum += var[i];

    return mean_sum / n;
}

// ------------------------------------------//
// ------------------------------------------//

double calculate_covariance(double *var_k, double *var_l, long int n, double mean_k, double mean_l)
{
    int chunk = CHUNKSIZE;

    register long int i;

    double cov_sum = 0.;

#pragma omp for schedule(dynamic, chunk)
    for (i = 0; i < n; i++)
        cov_sum += (var_k[i] - mean_k) * (var_l[i] - mean_l);

    return cov_sum / (n - 1);
}
