#include <stdio.h>
#include <stdlib.h>

#define N_SAMPLES	10
#define N_COEFFS	3

double	sample[N_SAMPLES] = {1,2,3,4,5,6,7,8,9,10};
double	coeff[N_COEFFS]= {0.5, 1, 0.5};
double	result[N_SAMPLES];

int data[5][5] =
{
    { 0, 0, 0, 0, 0 },
    { 0, 0, 1, 0, 0 },
    { 0, 0, 1, 0, 0 },
    { 0, 0, 1, 0, 0 },
    { 0, 0, 0, 0, 0 },
};

#define SIZEOF(a) (sizeof(a) / sizeof(a[0]))
int** Construct_New_Data(int (*old_data)[5], int rows, int cols){
    int** new_data = NULL;

    new_data = (int **)malloc(sizeof(int *)*(rows+2));
    for (int i=0; i<(rows+2); i++){
        new_data[i] = (int *)malloc(sizeof(int)*(cols+2));
    }

    
    for (int row = 0; row < (rows+2); row++) {
        for (int col = 0; col < (cols+2); col++) {
            if ( row == 0 || row == (rows+1)) {
                new_data[row][col] = 0;
            }
            else{
                if (col == 0 || col == (cols+1)) {
                    new_data[row][col] = 0;
                }
                else {
                    new_data[row][col] = old_data[row-1][col-1];
                }
            }

            printf("row:%d, col:%d value:%d\n",row, col, new_data[row][col] );
        }
    }

	int newrows = SIZEOF(new_data);
	int newcols = SIZEOF(new_data[0]);
	printf("row:%d, col:%d \n",newrows, newcols);
    for (int i = 0; i < rows+2; ++i)
	{
		for (int j = 0; j < cols+2; ++j)
		{
			printf("%d, ", new_data[i][j]);
		}
		printf("\n");
	}
    return new_data;
}


int main(int argc, char *arvg[])
{
	int old_row = SIZEOF(data);
    int old_col = SIZEOF(data[0]);
	int** new_data = Construct_New_Data(data,old_row,old_col);
	printf("%d\n", sizeof(int) );
	sizeof(char);
	return 0;
}
