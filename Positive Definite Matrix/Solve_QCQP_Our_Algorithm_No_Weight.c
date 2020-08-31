#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Problem Size:
#define N 1024 //number of decision variables
#define M 1 //number of constraints

//Algorithm Parameters:
#define MAX_ITER 10000000
#define TOLERANCE 0.001
#define LARGE 100






int main (int argc, char *argv[]){
    time_t start, end;
    int i, m, n_row, n_col;
    int nzcnt, k;
    FILE *fp1, *fp2, *fp3;
    char filename1[100];
    char filename2[100];
    char filename3[100];

    //Initiate Decision Variables:
    //Primal Variable:
    double *y = (double *)calloc(N, sizeof(double));
    double *x = (double *)calloc(N, sizeof(double));
    for (n_col = 0; n_col < N; n_col++){
        x[n_col] = 0.0;
    }
    //Dual Variable:
    double *mu = (double *)calloc(M, sizeof(double));
    double *lam = (double *)calloc(M, sizeof(double));
    for (m = 0; m < M; m++){
        lam[m] = 1.0;
    }

    //Read in Problem Formulation:
    double *q_obj = (double *)calloc(N, sizeof(double));
    double **q = (double **)calloc(M, sizeof(double *));
    for (m = 0; m < M; m++){
        q[m] = (double *)calloc(N, sizeof(double));
    }
    double **P_obj = (double **)calloc(N, sizeof(double *));
    for (n_row = 0; n_row < N; n_row++){
        P_obj[n_row] = (double *)calloc(N, sizeof(double));
    }
    double ***P = (double ***)calloc(M, sizeof(double **));
    for (m = 0; m < M; m++){
        P[m] = (double **)calloc(N, sizeof(double *));
        for (n_row = 0; n_row < N; n_row++){
            P[m][n_row] = (double *)calloc(N, sizeof(double));
        }
    }
    double *r = (double *)calloc(M, sizeof(double));
    
    //Initiate Objecive Function:
    //Linear Part:
    sprintf(filename1, "%d_%d/q_obj_%d_%d", M, N, M, N); //q_obj_M_N is the file containing the linear part of objective function
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < N; n_col++){
        fscanf(fp1, "%lf", &q_obj[n_col]);
    }
    fclose(fp1);
    //Quadratic Part:
    sprintf(filename1, "%d_%d/P_obj_%d_%d", M, N, M, N); //P_obj_M_N is the file containing the quadratic part of objective function
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    fscanf(fp1, "%d", &nzcnt);
    for (n_col = 0; n_col < N; n_col++){
        for (n_row = 0; n_row < N; n_row++){
            fscanf(fp1, "%lf", &P_obj[n_row][n_col]);
        }
    }
    fclose(fp1);
    //Initiate Constraints:
    for (m = 0; m < M; m++){
        //Linear Part:
        sprintf(filename1, "%d_%d/q_%d_%d_%d", M, N, m, M, N); //q_m_M_N is the file containing the linear part of constraint m
        if(!(fp1 = fopen(filename1,"r"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        fscanf(fp1, "%d", &nzcnt);
        for (n_col = 0; n_col < N; n_col++){
            fscanf(fp1, "%lf", &q[m][n_col]);
        }
        fclose(fp1);
        //Quadratic Part:
        sprintf(filename1, "%d_%d/P_%d_%d_%d", M, N, m, M, N); //P_m_M_N is the file containing the quadratic part of constraint m
        if(!(fp1 = fopen(filename1,"r"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        fscanf(fp1, "%d", &nzcnt);
        for (n_col = 0; n_col < N; n_col++){
            for (n_row = 0; n_row < N; n_row++){
                fscanf(fp1, "%lf", &P[m][n_row][n_col]);
            }
        }
        fclose(fp1);
    }
    //Right-Hand Side:
    sprintf(filename1, "%d_%d/r_%d_%d", M, N, M, N); //r_M_N is the file containing all the scalars of m constraints
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (m = 0; m < M; m++){
        fscanf(fp1, "%lf", &r[m]);
    }
    fclose(fp1);
printf("Initiation of Objective Function and Constraints is OK!\n");

    
    //Calculate Matrix Norm:
    double norm_P_obj = 0.0;
    double *norm_P = (double *)calloc(M, sizeof(double));
    double norm_PP = 0.0;
    for (n_row = 0; n_row < N; n_row++){
        for (n_col = 0; n_col < N; n_col++){
            norm_P_obj += pow(P_obj[n_row][n_col], (double)2);
        }
    }
    norm_P_obj = sqrt(norm_P_obj);
    for (m = 0; m < M; m++){
        norm_P[m] = 0.0;
        for (n_row = 0; n_row < N; n_row++){
            for (n_col = 0; n_col < N; n_col++){
                norm_P[m] += pow(P[m][n_row][n_col], (double)2);
            }
        }
        norm_PP += norm_P[m];
        norm_P[m] = sqrt(norm_P[m]);
    }
    norm_PP = sqrt(norm_PP);
    double norm_Q = 0.0;
    for (m = 0; m < M; m++){
        for (n_col = 0; n_col < N; n_col++){
            norm_Q += pow(q[m][n_col], (double)2);
        }
    }
    norm_Q = sqrt(norm_Q);
    double norm_x;
printf("Calculation of Matrix Norm is OK!\n");




    //Our Algorithm:
    double *rho = (double *)calloc(6, sizeof(double));
    double rho2_tmp;
    double *grad = (double *)calloc(N, sizeof(double));
    double *res = (double *)calloc(M, sizeof(double));
    double *epsilon = (double *)calloc(6, sizeof(double));
    double *weight = (double *)calloc(6, sizeof(double));
    for (i = 1; i < 6; i++){
        weight[i] = 1.0;
    }
    double *tol = (double *)calloc(2, sizeof(double));

    double a, b, c;

    sprintf(filename1, "%d_%d/OUT_Solve_QCQP_Our_Algorithm_No_Weight_Tolerance_%d_%d", M, N, M, N); //OUT_Solve_QCQP_Our_Algorithm_No_Weight_Tolerance_M_N is the file containing the tolerance of each iteration
    if(!(fp1 = fopen(filename1,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    sprintf(filename2, "%d_%d/OUT_Solve_QCQP_Our_Algorithm_No_Weight_Epsilon_%d_%d", M, N, M, N); //OUT_Solve_QCQP_Our_Algorithm_No_Weight_Epsilon_M_N is the file containing the epsilons of each iteration
    if(!(fp2 = fopen(filename2,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    sprintf(filename3, "%d_%d/OUT_Solve_QCQP_Our_Algorithm_No_Weight_Rho_%d_%d", M, N, M, N); //OUT_Solve_QCQP_Our_Algorithm_No_Weight_Rho_M_N is the file containing the rhos of each iteration
    if(!(fp3 = fopen(filename3,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }

    //Begin Iterations:
    start = time(NULL);
    for (k = 0; k < MAX_ITER; k++){
printf("This is the %d th iteration:\t", k);

        //Update epsilon:
        //epsilon[0] = 1/sqrt((double)(k+2));
        epsilon[0] = 0.0;
        for (i = 1; i < 6; i++){
            epsilon[i] = (1 - epsilon[0])/5;
        }
        fprintf(fp2, "%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n", epsilon[1], epsilon[2], epsilon[3], epsilon[4], epsilon[5]);

        //Update Gradient:
        tol[0] = 0.0;
        for (n_col = 0; n_col < N; n_col++){
            grad[n_col] = 0.0; 
            for (n_row = 0; n_row < N; n_row++){
                grad[n_col] += P_obj[n_col][n_row]*x[n_row];
            }
            grad[n_col] += q_obj[n_col];           
            for (m = 0; m < M; m++){
                for (n_row = 0; n_row < N; n_row++){
                    grad[n_col] += lam[m]*P[m][n_col][n_row]*x[n_row];
                }
                grad[n_col] += lam[m]*q[m][n_col];
            }
            tol[0] += pow(grad[n_col], (double)2);  
        }
        a = sqrt(tol[0]);
        tol[0] = sqrt(tol[0]/N);

        //Update Residual:
        tol[1] = 0.0;
        for (m = 0; m < M; m++){
            res[m] = 0.0;
            for (n_col = 0; n_col < N; n_col++){
                for (n_row = 0; n_row < N; n_row++){
                    res[m] += 0.5*x[n_col]*P[m][n_col][n_row]*x[n_row];
                }
                res[m] += q[m][n_col]*x[n_col];
            }
            res[m] += r[m];
            if (res[m] >= 0){
                tol[1] += pow(lam[m]*res[m], (double)2);
            }
            else{
                tol[1] += pow(lam[m]*(-res[m]), (double)2);
            }
        }
        tol[1] = sqrt(tol[1]/M);

        //Calculate Tolerance:
        fprintf(fp1, "%.16lf\t%.16lf\n", tol[0], tol[1]);
printf("%.16lf\t%.16lf\n", tol[0], tol[1]);
        if ((tol[0] < TOLERANCE) && (tol[1] < TOLERANCE)){
            break;
        }

        //Calculate Norm of x:
        norm_x = 0.0;
        for (n_col = 0; n_col < N; n_col++){
            norm_x += pow(x[n_col], (double)2);
        }
        norm_x = sqrt(norm_x);

        //Update Stepsize:
        //Rule 3:
        b = 2*norm_x;
        c = 2*epsilon[3]/norm_PP;
        if (a > TOLERANCE){
            rho[3] = (sqrt(pow(b, (double)2) + 4*a*c) - b)/(2*a);
        }
        else{
            if (b == 0.0){
                rho[3] = 2*epsilon[3];
            }
            else{
                rho[3] = c/b;
            }
        }
        //Rule 1:
        if (norm_P_obj == 0.0){
            rho[1] = epsilon[1];
        }
        else{
            rho[1] = epsilon[1]/norm_P_obj;
        }
        //Rule 2:
        rho[2] = LARGE;
        for (m = 0; m < M; m++){
            if (res[m] < 0){
                a = -res[m];
            }
            else{
                a = res[m];
            }
            b = lam[m];
            if (norm_P[m] == 0.0){
                c = epsilon[2]/M;
            }
            else{
                c = epsilon[2]/(M*norm_P[m]);
            }

            if (a > TOLERANCE){
                rho2_tmp = (sqrt(pow(b, (double)2) + 4*a*c) - b)/(2*a);
            }
            else{
                if (b == 0.0){
                    rho2_tmp = LARGE;
                }
                else{
                    rho2_tmp = c/b;
                }
            }

            if (rho2_tmp < rho[2]){
                rho[2] = rho2_tmp;
            }
        }
        //Rule 4:
        if (norm_Q == 0.0){
            rho[4] = epsilon[4];
        }
        else{
            rho[4] = epsilon[4]/norm_Q;
        }
        //Rule 5:
        if (norm_x == 0.0){
            rho[5] = epsilon[5];
        }
        else{
            rho[5] = epsilon[5]/(norm_x*norm_PP);
        }
        //Compare:
        rho[0] = rho[1];
        if (rho[2] < rho[0]){
            rho[0] = rho[2];
        }
        if (rho[3] < rho[0]){
            rho[0] = rho[3];
        }
        if (rho[4] < rho[0]){
            rho[0] = rho[4];
        }
        if (rho[5] < rho[0]){
            rho[0] = rho[5];
        }
        fprintf(fp3, "%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n", rho[1], rho[2], rho[3], rho[4], rho[5]);
//printf("%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n", rho[1], rho[2], rho[3], rho[4], rho[5]);
        //Adjust Weight:
        for (i = 1; i < 6; i++){
            weight[i] = weight[i]*(rho[0]/rho[i]);
        }

        //Update Primal and Dual Predictor:
        //Primal
        for (n_col = 0; n_col < N; n_col++){
            y[n_col] = x[n_col] - rho[0]*grad[n_col];
        }
        //Dual
        for (m = 0; m < M; m++){
            mu[m] = lam[m] + rho[0]*res[m];
            if (mu[m] < 0.0){
                mu[m] = 0.0;
            }
        }

        //Update Gradient:
        for (n_col = 0; n_col < N; n_col++){
            grad[n_col] = 0.0;
            for (n_row = 0; n_row < N; n_row++){
                grad[n_col] += P_obj[n_col][n_row]*y[n_row];
            }
            grad[n_col] += q_obj[n_col];
            for (m = 0; m < M; m++){
                for (n_row = 0; n_row < N; n_row++){
                    grad[n_col] += mu[m]*P[m][n_col][n_row]*y[n_row];
                }
                grad[n_col] += mu[m]*q[m][n_col];
            }
        }

        //Update Residual:
        for (m = 0; m < M; m++){
            res[m] = 0.0;
            for (n_col = 0; n_col < N; n_col++){
                for (n_row = 0; n_row < N; n_row++){
                    res[m] += 0.5*y[n_col]*P[m][n_col][n_row]*y[n_row];
                }
                res[m] += q[m][n_col]*y[n_col];
            }
            res[m] += r[m];
        }

        //Update Primal and Dual Corrector:
        //Primal
        for (n_col = 0; n_col < N; n_col++){
            x[n_col] += -rho[0]*grad[n_col];
        }
        //Dual
        for (m = 0; m < M; m++){
            lam[m] += rho[0]*res[m];
            if (lam[m] < 0.0){
                lam[m] = 0.0;
            }
        }
    }
    end = time(NULL);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    printf("Time use is %.16lf.\n", (double)(end-start));


    //Calculate Objective Function Value:
    double objval = 0.0;
    for (n_col = 0; n_col < N; n_col++){
        for (n_row = 0; n_row < N; n_row++){
            objval += 0.5*x[n_col]*P_obj[n_col][n_row]*x[n_row];
        }
        objval += q_obj[n_col]*x[n_col];
    }
    //Output Solution:
    sprintf(filename1, "%d_%d/OUT_Solve_QCQP_Our_Algorithm_No_Weight_%d_%d", M, N, M, N); //OUT_Solve_QCQP_Our_Algorithm_No_Weight_M_N is the file containing the solution
    if(!(fp1 = fopen(filename1,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    fprintf(fp1, "Time use is %.16lf.\n", (double)(end-start));
    fprintf(fp1, "The Objective Function Value is %.16lf.\n", objval);
    fprintf(fp1, "The Solution is:\n");
    for (n_col = 0; n_col < N; n_col++){
        fprintf(fp1, "%.16lf\n", x[n_col]);
    }
    fclose(fp1);
}
