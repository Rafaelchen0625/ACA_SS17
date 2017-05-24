#include <stdio.h>

#define N_SAMPLES	5
#define N_COEFFS	3

double	sample[N_SAMPLES] = {1, 2, 1, 2, 1};
double	coeff[N_COEFFS]= {0.5, 1, 0.5};
double	result[N_SAMPLES];

void smooth(double sample[], double coeff[], double result[], int n)
{
	int i, j;
	double norm=0.0;

	for (i=0; i<n; i++){
		if (i==0 || i==n-1){
			result[i] = sample[i];
		}else{
			result[i]=0.0;
			for (j=0; j<N_COEFFS; j++){
				if ( i == 1 ){
					norm+= coeff[j]>0 ? coeff[j] : -coeff[j];
				}
				result[i] += sample[i-1+j]*coeff[j];
			}
			result[i]/=norm;
		}
		printf("%f\n", result[i]);
	}
}

int main(int argc, char *arvg[])
{
	int i;

	if (N_SAMPLES>=N_COEFFS)
		smooth(sample, coeff, result, N_SAMPLES);
}
