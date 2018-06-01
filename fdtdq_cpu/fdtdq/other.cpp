#include "stdafx.h"
#include "other.h"
double *** generateZeroData(int nx, int ny, int nz)
{
	int i, j, k;
	double *** data = (double ***)malloc((nx) * sizeof(double **));
	for (i = 0; i<(nx); i++)
	{
		data[i] = (double **)malloc((ny) * sizeof(double *));
		for (j = 0; j<(ny); j++)
		{
			data[i][j] = (double *)malloc((nz) * sizeof(double));
			for (k = 0; k<(nz); k++)
			{
				data[i][j][k] = 0.0;
			}
		}
	}
	return data;
}
int setnsize(int*n, double width, double lambda_25)
{
	double dn = width / *n;
	*n = 3;
	while (dn >= lambda_25)
	{
		(*n)++;
		dn = width / (*n);
	}
	return 0;
}
double electric(int N)
{
	double dr = EARTH_HD / RE;
	double dt = dr / (2 * LIGHT_SPEED);
	double rtau = 0.3 / FREQUENCY;
	double tau = rtau / dt;
	return (double)(1000 * exp(-(N - 3 * tau)*(N - 3 * tau) / tau / tau));
}
