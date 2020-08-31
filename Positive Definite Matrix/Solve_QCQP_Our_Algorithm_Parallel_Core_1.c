#include <mpi.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Problem Size:
#define N 16384 //number of decision variables
#define M 1 //number of constraints
#define CONDITION_NUMBER 1.25
#define T1 1//number of threads for primal variables
#define T2 1//number of threads for dual variables

//Algorithm Parameters:
#define MAX_ITER 100000000
#define TOLERANCE 0.000001
#define LARGE 100.0






int main (int argc, char *argv[]){
    int size, rank;
    MPI_Comm comm_1, comm_2;
    MPI_Status status;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);






    //Declare All Variables:
    int size_1, rank_1;
    int size_2, rank_2;
    int n_t = N/T1; //number of primal variables per thread
    int m_t = M/T2; //number of dual variables per thread
    int t1;
    int t2;
    int i, m, n_row, n_col, k;
    double *y;
    double *x;
    double *mu;
    double *lam;
    double *q_obj;
    double **q;
    double **P_obj;
    double ***P;
    double *r;
    double *temp;
    FILE *fp1;
    char filename1[100];
    int nzcnt;
    double norm_Q_part;
    double norm_Q;
    double norm_P_obj_part;
    double norm_P_obj;
    double *norm_P_part;
    double *norm_P;
    double norm_PP;
    double norm_x_part;
    double norm_x;
    double *grad;    
    double *matrix_primal_part_1;
    double *matrix_primal;
    double *matrix_primal_part_2;
    double *res;
    double dual;
    double primal_matrix_primal_part;
    double primal_matrix_primal;
    double vector_primal_part;
    double vector_primal;
    double complm_part; 
    double complm;   
    double sq_grad_part;
    int flag = 888;      
    double *epsilon = (double *)calloc(6, sizeof(double));
    double *weight;
    double a, b, c;
    double rho_2_tmp, rho_2_min;
    double *rho = (double *)calloc(6, sizeof(double));            
    double *tol;
    double elapsed_time;
    double objval_part;
    double objval;
    double *solution;



    


    //Initiate Primal Communicator
    MPI_Comm_split(MPI_COMM_WORLD, !(rank < T1) ? MPI_UNDEFINED : 0, rank, &comm_1);
    if (rank < T1){
        MPI_Comm_size(comm_1, &size_1);
        MPI_Comm_rank(comm_1, &rank_1);
//printf("WORLD RANK/SIZE: %d/%d \t COMM1 RANK/SIZE: %d/%d\n", rank, size, rank_1, size_1);

        //Initiate Primal Decision Variables:
        y = (double *)calloc(n_t, sizeof(double));
        x = (double *)calloc(n_t, sizeof(double));
        for (n_col = 0; n_col < n_t; n_col++){
            x[n_col] = 0.0;
        }

        //Read in Problem Formulation:
        q_obj = (double *)calloc(n_t, sizeof(double));
        q = (double **)calloc(M, sizeof(double *));
        for (m = 0; m < M; m++){
            q[m] = (double *)calloc(n_t, sizeof(double));
        }
        P_obj = (double **)calloc(N, sizeof(double *));
        for (n_row = 0; n_row < N; n_row++){
            P_obj[n_row] = (double *)calloc(n_t, sizeof(double));
        }
        P = (double ***)calloc(M, sizeof(double **));
        for (m = 0; m < M; m++){
            P[m] = (double **)calloc(N, sizeof(double *));
            for (n_row = 0; n_row < N; n_row++){
                P[m][n_row] = (double *)calloc(n_t, sizeof(double));
            }
        }

        //Initiate Objective Function:
        //Linear Part:       
        if (!(rank_1)){                                    
            sprintf(filename1, "%d_%d/Condition_Number_%.2f/q_obj_%d_%d", M, N, CONDITION_NUMBER, M, N); //q_obj_M_N is the file containing the linear objective function
            if(!(fp1 = fopen(filename1,"r"))){
                printf("Cannot open this file\n");
                exit(0);
            }
            temp = (double *)calloc(N, sizeof(double));
            for (n_col = 0; n_col < N; n_col++){
                fscanf(fp1, "%lf", &temp[n_col]);
            }
            fclose(fp1);
        }
        MPI_Scatter(temp, n_t, MPI_DOUBLE, q_obj, n_t, MPI_DOUBLE, 0, comm_1);
        //Quadratic Part:
        if (!(rank_1)){           
            sprintf(filename1, "%d_%d/Condition_Number_%.2f/P_obj_%d_%d", M, N, CONDITION_NUMBER, M, N); //P_obj_M_N is the file containing the quadratic objective function
            if(!(fp1 = fopen(filename1,"r"))){
                printf("Cannot open this file\n");
                exit(0);
            }
            fscanf(fp1, "%d", &nzcnt);
        }
        for (n_row = 0; n_row < N; n_row++){
            if (!(rank_1)){
                for (n_col = 0; n_col < N; n_col++){
                    fscanf(fp1, "%lf", &temp[n_col]);
                }
            }
            MPI_Scatter(temp, n_t, MPI_DOUBLE, P_obj[n_row], n_t, MPI_DOUBLE, 0, comm_1);
        }
        if (!(rank_1)){
            fclose(fp1);
        }
        //Initiate Constraints:
        for (m = 0; m < M; m++){
            //Linear Part:
            if (!(rank_1)){
                sprintf(filename1, "%d_%d/Condition_Number_%.2f/q_%d_%d_%d", M, N, CONDITION_NUMBER, m, M, N); //q_m_M_N is the file containing the linear part of constraint m
                if(!(fp1 = fopen(filename1,"r"))){
                    printf("Cannot open this file\n");
                    exit(0);
                }
                fscanf(fp1, "%d", &nzcnt);
                for (n_col = 0; n_col < N; n_col++){
                    fscanf(fp1, "%lf", &temp[n_col]);
                }
                fclose(fp1);
            }
            MPI_Scatter(temp, n_t, MPI_DOUBLE, q[m], n_t, MPI_DOUBLE, 0, comm_1);
            //Quadratic Part:
            if (!(rank_1)){                
                sprintf(filename1, "%d_%d/Condition_Number_%.2f/P_%d_%d_%d", M, N, CONDITION_NUMBER, m, M, N); //P_m_M_N is the file containing the quadratic part of constraint m
                if(!(fp1 = fopen(filename1,"r"))){
                    printf("Cannot open this file\n");
                    exit(0);
                }
                fscanf(fp1, "%d", &nzcnt);
            }
            for (n_row = 0; n_row < N; n_row++){
                if (!(rank_1)){
                    for (n_col = 0; n_col < N; n_col++){
                        fscanf(fp1, "%lf", &temp[n_col]);
                    }
                }
                MPI_Scatter(temp, n_t, MPI_DOUBLE, P[m][n_row], n_t, MPI_DOUBLE, 0, comm_1);
            }
            if (!(rank_1)){
                fclose(fp1);
            }
        }       
        if (!(rank_1)){
        free(temp);
printf("Initiation of Objective Function and Left-Hand Side of Constraints is OK!\n");
        }

        //Calculate Matrix Norm:
        norm_Q_part = 0.0;
        for (m = 0; m < M; m++){
            for (n_col = 0; n_col < n_t; n_col++){
                norm_Q_part += pow(q[m][n_col], (double)2);
            }
        }
        MPI_Reduce(&norm_Q_part, &norm_Q, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);  
        norm_P_obj_part = 0.0;
        for (n_row = 0; n_row < N; n_row++){
            for (n_col = 0; n_col < n_t; n_col++){
                norm_P_obj_part += pow(P_obj[n_row][n_col], (double)2);
            }
        }
        MPI_Reduce(&norm_P_obj_part, &norm_P_obj, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        if (!(rank_1)){
            norm_P = (double *)calloc(M, sizeof(double));
        }
        norm_P_part = (double *)calloc(M, sizeof(double));
        for (m = 0; m < M; m++){
            norm_P_part[m] = 0.0;
            for (n_row = 0; n_row < N; n_row++){
                for (n_col = 0; n_col < n_t; n_col++){
                    norm_P_part[m] += pow(P[m][n_row][n_col], (double)2);
                }
            }
            MPI_Reduce(&norm_P_part[m], &norm_P[m], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        }
        if (!(rank_1)){
            norm_Q = sqrt(norm_Q);
            norm_P_obj = sqrt(norm_P_obj);
            norm_PP = 0.0;
            for (m = 0; m < M; m++){
                norm_PP += norm_P[m];
                norm_P[m] = sqrt(norm_P[m]);
            }
            norm_PP = sqrt(norm_PP);
//printf("%.16lf\t%.16lf\t%.16lf\n", norm_Q, norm_P_obj, norm_PP);
printf("Calculation of Matrix Norms is OK!\n");
        }

        //Initiate Our Algorithm Variables:
        grad = (double *)calloc(n_t, sizeof(double));
        matrix_primal_part_1 = (double *)calloc(N, sizeof(double));
        matrix_primal_part_2 = (double *)calloc(n_t, sizeof(double));
        if (!(rank_1)){
        	matrix_primal = (double *)calloc(N, sizeof(double));
            weight = (double *)calloc(6, sizeof(double));
            for (i = 1; i < 6; i++){
                weight[i] = 1.0;
            }        
            tol = (double *)calloc(2, sizeof(double));

            //Open File
            sprintf(filename1, "%d_%d/Condition_Number_%.2f/OUT_Solve_QCQP_Our_Algorithm_Parallel_Tolerance_%d_%d_Core_%d", M, N, CONDITION_NUMBER, M, N, T1+T2); //OUT_Solve_QCQP_Our_Algorithm_Parallel_Tolerance_M_N_T1+T2 is the file containing the tolerance of each iteration
            if(!(fp1 = fopen(filename1,"w"))){
                printf("Cannot open this file\n");
                exit(0);
            }       
        }
    }


    //Initiate Dual Communicator
    MPI_Comm_split(MPI_COMM_WORLD, (rank < T1) ? MPI_UNDEFINED : 0, rank, &comm_2);
    if (rank >= T1){
        MPI_Comm_size(comm_2, &size_2);
        MPI_Comm_rank(comm_2, &rank_2);
//printf("WORLD RANK/SIZE: %d/%d \t COMM2 RANK/SIZE: %d/%d\n", rank, size, rank_2, size_2);

        //Initiate Dual Decision Variables:
        mu = (double *)calloc(m_t, sizeof(double));
        lam = (double *)calloc(m_t, sizeof(double));
        for (m = 0; m < m_t; m++){
            lam[m] = 1.0;
        }

        //Read in Problem Formulation:
        r = (double *)calloc(m_t, sizeof(double));

        //Initiate Constraints:
        //Right-Hand Side:
        if (!(rank_2)){            
            sprintf(filename1, "%d_%d/Condition_Number_%.2f/r_%d_%d", M, N, CONDITION_NUMBER, M, N); //r_M_N is the file containing all the scalars of m constraints
            if(!(fp1 = fopen(filename1,"r"))){
                printf("Cannot open this file\n");
                exit(0);
            }
            temp = (double *)calloc(M, sizeof(double));
            for (m = 0; m < M; m++){
                fscanf(fp1, "%lf", &temp[m]);
            }
            fclose(fp1);
        }
        MPI_Scatter(temp, m_t, MPI_DOUBLE, r, m_t, MPI_DOUBLE, 0, comm_2);
        if (!(rank_2)){            
            free(temp);
printf("Initiation of Right-Hand Side of Constraints is OK!\n");
        }

        //Initiate Matrix Norm:
        norm_P_part = (double *)calloc(m_t, sizeof(double));
        if (!(rank_2)){
        	norm_P = (double *)calloc(M, sizeof(double));
        }

        //Initiate Our Algorithm Variables:
        res = (double *)calloc(m_t, sizeof(double)); 
    }


    //Copy Matrix Norm from Primal to Dual Communicators:
    MPI_Barrier(MPI_COMM_WORLD);
    //Send Matrix Norm:
    if (!(rank)){
    	MPI_Send(norm_P, M, MPI_DOUBLE, T1, 88, MPI_COMM_WORLD);
    }
    //Receive Matrix Norm:
    if (rank == T1){
    	MPI_Recv(norm_P, M, MPI_DOUBLE, 0, 88, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    //Scatter Matrix Norm:
    if (rank >= T1){
    	MPI_Scatter(norm_P, m_t, MPI_DOUBLE, norm_P_part, m_t, MPI_DOUBLE, 0, comm_2);
    	if (!(rank_2)){
    		free(norm_P);
    	}
    }



    


    //Begin Iterations:
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    for (k = 0; k < MAX_ITER; k++){
/*
if (!(rank)){
    printf("This is the %d th iteration:\t", k);
}
*/
        //Update Gradient and Residual:
        if (rank < T1){
            //Calculate P_obj*x:
        	for (n_row = 0; n_row < N; n_row++){
        		matrix_primal_part_1[n_row] = 0.0;
        		for (n_col = 0; n_col < n_t; n_col++){
        			matrix_primal_part_1[n_row] += P_obj[n_row][n_col]*x[n_col];
        		}
        		MPI_Reduce(&matrix_primal_part_1[n_row], &matrix_primal[n_row], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	}
        	MPI_Scatter(matrix_primal, n_t, MPI_DOUBLE, matrix_primal_part_2, n_t, MPI_DOUBLE, 0, comm_1);
        	//Update Gradient with P_obj*x + q_obj:
        	for (n_col = 0; n_col < n_t; n_col++){
        		grad[n_col] =  matrix_primal_part_2[n_col] + q_obj[n_col];
        	}
        }
        else{
            complm_part = 0.0;
        }
        for (m = 0; m < M; m++){
            //Send Dual:
        	if (rank == T1 + (m/m_t)){
        		MPI_Send(&lam[(m % m_t)], 1, MPI_DOUBLE, 0, 0*M + m, MPI_COMM_WORLD);
        	}
        	//Receive Dual:
        	if (!(rank)){
        		MPI_Recv(&dual, 1, MPI_DOUBLE, (T1 + m/m_t), 0*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	}
        	if (rank < T1){
        		//Broadcast Dual:
        		MPI_Bcast(&dual, 1, MPI_DOUBLE, 0, comm_1);
        		//Calculate P_m*x:
        		for (n_row = 0; n_row < N; n_row++){
        			matrix_primal_part_1[n_row] = 0.0;
        			for (n_col = 0; n_col < n_t; n_col++){
        				matrix_primal_part_1[n_row] += P[m][n_row][n_col]*x[n_col];
        			}
        			MPI_Reduce(&matrix_primal_part_1[n_row], &matrix_primal[n_row], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        		}
        		MPI_Scatter(matrix_primal, n_t, MPI_DOUBLE, matrix_primal_part_2, n_t, MPI_DOUBLE, 0, comm_1);
        		//Update Gradient with lam_m*(P_m*x + q_m):
        		for (n_col = 0; n_col < n_t; n_col++){
        			grad[n_col] += dual*(matrix_primal_part_2[n_col] + q[m][n_col]);
        		}
        		//Calculate x*P_m*x and q_m*x:
        		primal_matrix_primal_part = 0.0;
        		vector_primal_part = 0.0;
        		for (n_col = 0; n_col < n_t; n_col++){
        			primal_matrix_primal_part += matrix_primal_part_2[n_col]*x[n_col];
        			vector_primal_part += q[m][n_col]*x[n_col];
        		}
        		MPI_Reduce(&primal_matrix_primal_part, &primal_matrix_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        		MPI_Reduce(&vector_primal_part, &vector_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	}
        	//Send 0.5*x*P_m*x + q_m*x:
        	if (!(rank)){
        		primal_matrix_primal = 0.5*primal_matrix_primal;
        		primal_matrix_primal += vector_primal;
        		MPI_Send(&primal_matrix_primal, 1, MPI_DOUBLE, (T1 + m/m_t), 1*M + m, MPI_COMM_WORLD);
        	}
        	//Receive 0.5*x*P_m*x + q_m*x and Update Residual with 0.5*x*P_m*x + q_m*x + r_m:
        	if (rank == T1 + (m/m_t)){
        		MPI_Recv(&res[(m % m_t)], 1, MPI_DOUBLE, 0, 1*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        		res[(m % m_t)] += r[(m % m_t)];
                if (res[(m % m_t)] >= 0){
                    complm_part += pow(lam[(m % m_t)]*res[(m % m_t)], (double)2);
                }
                else{
                    complm_part += pow(lam[(m % m_t)]*(-res[(m % m_t)]), (double)2);
                }
        	}
        }
        if (rank < T1){
            sq_grad_part = 0.0;
            for (n_col = 0; n_col < n_t; n_col++){
                sq_grad_part += pow(grad[n_col], (double)2);
            }
            MPI_Reduce(&sq_grad_part, &tol[0], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        }
        else{
            MPI_Reduce(&complm_part, &complm, 1, MPI_DOUBLE, MPI_SUM, 0, comm_2);
            if (rank == T1){
                MPI_Send(&complm, 1, MPI_DOUBLE, 0, 888, MPI_COMM_WORLD);
            }
        }
        if (!(rank)){
            MPI_Recv(&tol[1], 1, MPI_DOUBLE, T1, 888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            a = sqrt(tol[0]);
            
            //Calculate Tolerance:
            tol[0] = sqrt(tol[0]/N);
            tol[1] = sqrt(tol[1]/M);
            fprintf(fp1, "%.16lf\t%.16lf\n", tol[0], tol[1]);
//printf("%.16lf\t%.16lf\n", tol[0], tol[1]);
            if ((tol[0] < TOLERANCE) && (tol[1] < TOLERANCE)){
                flag = 0;
            }
        }
        MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!(flag)){
            break;
        }

        //Calculate Norm of x:
        if (rank < T1){
        	norm_x_part = 0.0;
        	for (n_col = 0; n_col < n_t; n_col++){
        		norm_x_part += pow(x[n_col], (double)2);
        	}
        	MPI_Reduce(&norm_x_part, &norm_x, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        }
        if (!(rank)){
        	norm_x = sqrt(norm_x);

            //Update epsilon:
            //epsilon[0] = 1/sqrt((double)(k+2));
            epsilon[0] = 0.0;
            weight[0] = 0.0;
            for (i = 1; i < 6; i++){
                weight[0] += weight[i];
            }
            for (i = 1; i < 6; i++){
                epsilon[i] = (1 - epsilon[0])*(weight[i]/weight[0]);
            }
            MPI_Send(&epsilon[2], 1, MPI_DOUBLE, T1, 8888, MPI_COMM_WORLD);
        }

        //Update Stepsize:
        if (rank >= T1){
            //Rule 2:
            if (rank == T1){
                MPI_Recv(&epsilon[2], 1, MPI_DOUBLE, 0, 8888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Bcast(&epsilon[2], 1, MPI_DOUBLE, 0, comm_2);
            rho[2] = LARGE;
            for (m = 0; m < m_t; m++){
                if (res[m] < 0){
                    a = -res[m];
                }
                else{
                    a = res[m];
                }
                b = lam[m];
                if (norm_P_part[m] == 0.0){
                    c = epsilon[2]/M;
                }
                else{
                    c = epsilon[2]/(M*norm_P_part[m]);
                }
                if (a > TOLERANCE){
                    rho_2_tmp = (sqrt(pow(b, (double)2) + 4*a*c) - b)/(2*a);
                }
                else{
                    if (b == 0.0){
                        rho_2_tmp = LARGE;
                    }
                    else{
                        rho_2_tmp = c/b;
                    }
                }
                if (rho_2_tmp < rho[2]){
                    rho[2] = rho_2_tmp;
                }
            }
            MPI_Reduce(&rho[2], &rho_2_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm_2);
            if (rank == T1){
                MPI_Send(&rho_2_min, 1, MPI_DOUBLE, 0, 88888, MPI_COMM_WORLD);
            }
        }
        if (!(rank)){
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
            MPI_Recv(&rho[2], 1, MPI_DOUBLE, T1, 88888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
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
            //Adjust Weight:
            for (i = 1; i < 6; i++){
                weight[i] = weight[i]*(rho[0]/rho[i]);
            }
//printf("%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n", rho[1], rho[2], rho[3], rho[4], rho[5]);
        }


        //Update Primal and Dual Predictor:
        MPI_Bcast(&rho[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //Primal
        if (rank < T1){
        	for (n_col = 0; n_col < n_t; n_col++){
        		y[n_col] = x[n_col] - rho[0]*grad[n_col];
        	}
        }       
        //Dual
        else{
        	for (m = 0; m < m_t; m++){
        		mu[m] = lam[m] + rho[0]*res[m];
        		if (mu[m] < 0.0){
        			mu[m] = 0.0;
        		}
        	}
        }


        //Update Gradient and Residual:
        if (rank < T1){
        	//Calculate P_obj*y:
        	for (n_row = 0; n_row < N; n_row++){
        		matrix_primal_part_1[n_row] = 0.0;
        		for (n_col = 0; n_col < n_t; n_col++){
        			matrix_primal_part_1[n_row] += P_obj[n_row][n_col]*y[n_col];
        		}
        		MPI_Reduce(&matrix_primal_part_1[n_row], &matrix_primal[n_row], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	}
        	MPI_Scatter(matrix_primal, n_t, MPI_DOUBLE, matrix_primal_part_2, n_t, MPI_DOUBLE, 0, comm_1);
        	//Update Gradient with P_obj*y + q_obj:
        	for (n_col = 0; n_col < n_t; n_col++){
        		grad[n_col] =  matrix_primal_part_2[n_col] + q_obj[n_col];
        	}
        }       
        for (m = 0; m < M; m++){
        	//Send Dual:
        	if (rank == T1 + (m/m_t)){
        		MPI_Send(&mu[(m % m_t)], 1, MPI_DOUBLE, 0, 2*M + m, MPI_COMM_WORLD);
        	}
        	//Receive Dual:
        	if (!(rank)){
        		MPI_Recv(&dual, 1, MPI_DOUBLE, (T1 + m/m_t), 2*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	}
        	if (rank < T1){
        		//Broadcast Dual:
        		MPI_Bcast(&dual, 1, MPI_DOUBLE, 0, comm_1);
        		//Calculate P_m*y:
        		for (n_row = 0; n_row < N; n_row++){
        			matrix_primal_part_1[n_row] = 0.0;
        			for (n_col = 0; n_col < n_t; n_col++){
        				matrix_primal_part_1[n_row] += P[m][n_row][n_col]*y[n_col];
        			}
        			MPI_Reduce(&matrix_primal_part_1[n_row], &matrix_primal[n_row], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        		}
        		MPI_Scatter(matrix_primal, n_t, MPI_DOUBLE, matrix_primal_part_2, n_t, MPI_DOUBLE, 0, comm_1);
        		//Update Gradient with mu_m*(P_m*y + q_m):
        		for (n_col = 0; n_col < n_t; n_col++){
        			grad[n_col] += dual*(matrix_primal_part_2[n_col] + q[m][n_col]);
        		}
        		//Calculate x*P_m*x and q_m*x:
        		primal_matrix_primal_part = 0.0;
        		vector_primal_part = 0.0;
        		for (n_col = 0; n_col < n_t; n_col++){
        			primal_matrix_primal_part += matrix_primal_part_2[n_col]*y[n_col];
        			vector_primal_part += q[m][n_col]*y[n_col];
        		}
        		MPI_Reduce(&primal_matrix_primal_part, &primal_matrix_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        		MPI_Reduce(&vector_primal_part, &vector_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	}
        	//Send 0.5*y*P_m*y + q_m*y:
        	if (!(rank)){
        		primal_matrix_primal = 0.5*primal_matrix_primal;
        		primal_matrix_primal += vector_primal;
        		MPI_Send(&primal_matrix_primal, 1, MPI_DOUBLE, (T1 + m/m_t), 3*M + m, MPI_COMM_WORLD);
        	}
        	//Receive 0.5*y*P_m*y + q_m*y and Update Residual with 0.5*y*P_m*y + q_m*y + r_m:
        	if (rank == T1 + (m/m_t)){
        		MPI_Recv(&res[(m % m_t)], 1, MPI_DOUBLE, 0, 3*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        		res[(m % m_t)] += r[(m % m_t)];
        	}
        }


        //Update Primal and Dual Corrector:
        //Primal
        if (rank < T1){
        	for (n_col = 0; n_col < n_t; n_col++){
        		x[n_col] += -rho[0]*grad[n_col];
        	}
        }       
        //Dual
        else{
        	for (m = 0; m < m_t; m++){
        		lam[m] += rho[0]*res[m];
        		if (lam[m] < 0.0){
        			lam[m] = 0.0;
        		}
        	}
        }      
    }
    elapsed_time += MPI_Wtime();
    if (!(rank)){
    	fclose(fp1);
    	printf("Time use is %.16lf.\n", elapsed_time);
    }

    //Calculate Objective Function Value:
    if (rank < T1){
    	for (n_row = 0; n_row < N; n_row++){
    		matrix_primal_part_1[n_row] = 0.0;
    		for (n_col = 0; n_col < n_t; n_col++){
    			matrix_primal_part_1[n_row] += P_obj[n_row][n_col]*x[n_col];
    		}
    		MPI_Reduce(&matrix_primal_part_1[n_row], &matrix_primal[n_row], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
    	}
    	MPI_Scatter(matrix_primal, n_t, MPI_DOUBLE, matrix_primal_part_2, n_t, MPI_DOUBLE, 0, comm_1);

    	objval_part = 0.0;
    	for (n_col = 0; n_col < n_t; n_col++){
    		objval_part += 0.5*matrix_primal_part_2[n_col]*x[n_col] + q_obj[n_col]*x[n_col];
    	}
    	if (!(rank_1)){
    		solution = (double *)calloc(N, sizeof(double));
    	}
    	MPI_Reduce(&objval_part, &objval, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
    	MPI_Gather(x, n_t, MPI_DOUBLE, solution, n_t, MPI_DOUBLE, 0, comm_1);

    	if (!(rank_1)){
    		//Output Solution:
    		sprintf(filename1, "%d_%d/Condition_Number_%.2f/OUT_Solve_QCQP_Our_Algorithm_Parallel_%d_%d_Core_%d", M, N, CONDITION_NUMBER, M, N, T1+T2); //OUT_Solve_QCQP_Our_Algorithm_Parallel_M_N_T1+T2 is the file containing the solution
    		if(!(fp1 = fopen(filename1,"w"))){
    			printf("Cannot open this file\n");
    			exit(0);
    		}
    		fprintf(fp1, "Time use is %.16lf.\n", elapsed_time);
    		fprintf(fp1, "The Objective Function Value is %.16lf.\n", objval);
    		fprintf(fp1, "The Solution is:\n");
    		for (n_col = 0; n_col < N; n_col++){
    			fprintf(fp1, "%.16lf\n", solution[n_col]);
    		}
    		fclose(fp1);
    	}
    }  
    //Terminate:
    if (rank < T1){
        MPI_Comm_free(&comm_1);
    }
    else{
        MPI_Comm_free(&comm_2);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
