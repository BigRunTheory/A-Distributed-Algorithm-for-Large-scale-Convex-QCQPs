#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Problem Size:
#define N 16384 //number of decision variables
#define M 1 //number of constraints
#define EIGENVALUE_L 4.0 //smallest eigenvale
#define EIGENVALUE_H 5.0 //largest eigenvalue
#define CONDITION_NUMBER 1.25 //condition number is EIGENVALUE_H/EIGENVALUE_L.
#define Q_S 1.0
#define R_S 1.0
//condition number: 1 (5.0, 5.0)
//condition number: 1.25 (4.0, 5.0)
//condition number: 10.00 (0.8, 8.0)
//condition number: 100.00 (0.1, 10.0)
//condition number: 10000.00 (0.003, 30.0)
//condition number: 1000000.00 (0.00002, 20.0)






int generate_random_vector (double *v, int dim, double lower, double upper){
    //Generate a vector of dim random doubles ranging from [lower, upper], stored in v.
    int i;
    int nz_cnt = 0;
    for (i = 0; i < dim; i++){
    	v[i] = (((double)rand())/RAND_MAX)*(upper - lower) + lower;
    	if (v[i] != 0.0){
    		nz_cnt++;
    	}
    }
    return (nz_cnt);
}






int generate_random_orthogonal_matrix (double **A, int dim){
	//Generate a dim by dim othorgonal matrix, stored in A.
	int i, j, k;
	int nz_cnt;
	double sum = 0.0;
	double *v = (double *)calloc(dim, sizeof(double));

	nz_cnt = generate_random_vector(v, dim, -1.0, 1.0);

	for (i = 0; i < dim; i++){
		sum += pow(v[i], (double)2);
	}
	sum = sqrt(sum);
	for (i = 0; i < dim; i++){
		v[i] = v[i]/sum;
	} 

	double cos_phi = v[0];
	sum = 0.0;
	for (i = 1; i < dim; i++){
		sum += pow(v[i], (double)2);
	}
	double sin_phi = sqrt(sum);
	for (i = 1; i < dim; i++){
		v[i] = v[i]/sin_phi;
	}

	for (i = 0; i < dim; i++){
		for (j = i; j < dim; j++){
			if (i == 0){
				if (j == i){
					A[i][j] = cos_phi;
				}
				else{
					A[i][j] = -v[j]*sin_phi;
					A[j][i] = v[j]*sin_phi;
				}
			}
			else{
				if (j == i){
					A[i][j] = 1 + pow(v[i], (double)2)*(cos_phi - 1);
				}
				else{
					A[i][j] = v[i]*v[j]*(cos_phi - 1);
					A[j][i] = A[i][j];
				}
			}
//printf("%.16lf\n", A[i][j]);
		}
	}
    
	int flag = 0;
	/*
	for (i = 0; i < dim; i++){
		for (j = 0; j < dim; j++){
			sum = 0.0;
			for (k = 0; k < dim; k++){
				sum += A[i][k]*A[j][k];
			}
			if (i == j){
				if (sum != 1.0){
					flag = 88;
				}
			}
			else{
				if (sum != 0.0){
					flag = 88;
				}
			}
printf("%.16lf\n", sum);
		}
	}
	*/

	return (flag);
}






int matrix_multiplication_diagonal_othorgonal (double *v, double **A, int dim, double **B){
	//Calculate B = A^T V A.
	int i, j, k;
	int nz_cnt = 0;

	for (i = 0; i < dim; i++){
		for (j = i; j < dim; j++){
			B[i][j] = 0.0;
			for (k = 0; k < dim; k++){
				B[i][j] += A[i][k]*v[k]*A[j][k];
			}
			if (B[i][j] != 0.0){
				nz_cnt++;
			}
			if (j != i){
				B[j][i] = B[i][j];
				if (B[j][i] != 0.0){
					nz_cnt++;
				}
			}
		}
	}

	return (nz_cnt);
}






int main (int argc, char *argv[]){
	srand(time(NULL));
	int m, n_row, n_col;
	int nz_cnt;
	int flag;


	double *q = (double *)calloc(N, sizeof(double));
	double *d = (double *)calloc(N, sizeof(double));
    double **P = (double **)calloc(N, sizeof(double *));
    for (n_row = 0; n_row < N; n_row++){
    	P[n_row] = (double *)calloc(N, sizeof(double));
    }
    double **O = (double **)calloc(N, sizeof(double *));
    for (n_row = 0; n_row < N; n_row++){
    	O[n_row] = (double *)calloc(N, sizeof(double));
    }

	FILE *fp1;
    char  filename1[100];


	//Create Objective Function:
	//Linear Part:
	nz_cnt = generate_random_vector (q, N, -Q_S, Q_S);
    sprintf(filename1, "%d_%d/Condition_Number_%.2f/q_obj_%d_%d", M, N, CONDITION_NUMBER, M, N); //q_obj_M_N is the file containing the linear objective function
    if(!(fp1 = fopen(filename1, "w"))){
		printf("Cannot open this file\n");
		exit(0);
	}
	for (n_col = 0; n_col < N; n_col++){
		fprintf(fp1, "%.16lf\n", q[n_col]);
	}
	fclose(fp1);

	//Quadratic Part:
	nz_cnt = generate_random_vector (d, N, EIGENVALUE_L, EIGENVALUE_H);
	d[0] = EIGENVALUE_L;
	d[N-1] = EIGENVALUE_H;
	flag = generate_random_orthogonal_matrix (O, N);
	if (flag != 0){
		printf("Error in generating othorgonal matrix!\n");
	}
	nz_cnt = matrix_multiplication_diagonal_othorgonal (d, O, N, P);
	sprintf(filename1, "%d_%d/Condition_Number_%.2f/P_obj_%d_%d", M, N, CONDITION_NUMBER, M, N); //P_obj_M_N is the file containing the quadratic objective function
    if(!(fp1 = fopen(filename1, "w"))){
		printf("Cannot open this file\n");
		exit(0);
	}
	fprintf(fp1, "%d\n", nz_cnt);
	for (n_col = 0; n_col < N; n_col++){
		for (n_row = 0; n_row < N; n_row++){
			fprintf(fp1, "%.16lf\n", P[n_row][n_col]);
		}
	}
	fclose(fp1);
printf("Creation of Objective Function is OK!\n");

    //Creat Constraints:
    for (m = 0; m < M; m++){
printf("%d\n", m);
    	//Linear Part:
    	nz_cnt = generate_random_vector (q, N, -Q_S, Q_S);
    	sprintf(filename1, "%d_%d/Condition_Number_%.2f/q_%d_%d_%d", M, N, CONDITION_NUMBER, m, M, N); //q_m_M_N is the file containing the linear part of constraint m
    	if(!(fp1 = fopen(filename1, "w"))){
    		printf("Cannot open this file\n");
    		exit(0);
    	}
    	fprintf(fp1, "%d\n", nz_cnt);
    	for (n_col = 0; n_col < N; n_col++){
    		fprintf(fp1, "%.16lf\n", q[n_col]);
    	}
    	fclose(fp1);

    	//Quadratic Part:
    	nz_cnt = generate_random_vector (d, N, EIGENVALUE_L, EIGENVALUE_H);
    	d[0] = EIGENVALUE_L;
    	d[N-1] = EIGENVALUE_H;
    	flag = generate_random_orthogonal_matrix (O, N);
    	if (flag){
    		printf("Error in generating othorgonal matrix!\n");
    	}
    	nz_cnt = matrix_multiplication_diagonal_othorgonal (d, O, N, P);
    	sprintf(filename1, "%d_%d/Condition_Number_%.2f/P_%d_%d_%d", M, N, CONDITION_NUMBER, m, M, N); //P_m_M_N is the file containing the quadratic part of constraint m
    	if(!(fp1 = fopen(filename1, "w"))){
    		printf("Cannot open this file\n");
    		exit(0);
    	}
    	fprintf(fp1, "%d\n", nz_cnt);
    	for (n_col = 0; n_col < N; n_col++){
    		for (n_row = 0; n_row < N; n_row++){
    			fprintf(fp1, "%.16lf\n", P[n_row][n_col]);
    		}
    	}
    	fclose(fp1);
    }

    //Right-Hand Side:
    nz_cnt = generate_random_vector (q, M, -R_S, 0.0);
    sprintf(filename1, "%d_%d/Condition_Number_%.2f/r_%d_%d", M, N, CONDITION_NUMBER, M, N); //r_M_N is the file containing all the scalars of m constraints
    if(!(fp1 = fopen(filename1,"w"))){
    	printf("Cannot open this file\n");
    	exit(0);
    }
    for (m = 0; m < M; m++){
    	fprintf(fp1, "%.16lf\n", q[m]);
    }
    fclose(fp1);
printf("Creation of %d Constraints is OK!\n", M);
}





