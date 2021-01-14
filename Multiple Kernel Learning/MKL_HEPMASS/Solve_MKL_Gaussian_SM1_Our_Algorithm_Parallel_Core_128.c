#include <mpi.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


//static char FOLDER[] = "2_norm";
static char FOLDER[] = "HEPMASS";
//static char FOLDER[] = "Ionosphere";
//Problem Size:
#define NUM_DATA 13760
#define N 11008 //number of decision variables
#define M 5 //number of constraints
#define T1 128//number of threads for primal variables
#define T2 5//number of threads for dual variables

//Algorithm Parameters:
#define MAX_ITER 110000
//#define MAX_ITER 100000000
#define TOLERANCE_1 0.001
#define TOLERANCE_2 0.001
#define LARGE 100.0
#define C 5.0






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
    int i, m, n_row, n_col, n_att, k;
    double *y;
    double *x;
    double *v;
    double *u;
    double *mu;
    double *lam;
    double *nu;
    double *gam;
    double *q_obj;
    double *c_obj;
    double **c_vec;
    double *temp;
    FILE *fp1;
    char filename1[100];
    double **A;
    int nzcnt;
    double ***P;        
    double *norm_P_part;
    double *norm_P;
    double norm_PP;
    double norm_C;
    double norm_A_part;
    double norm_A;
    double norm_x_part;
    double norm_x;
    double *grad_x;
    double *nc_grad_x;
    double *grad_u;
    double *res_ineq;
    double *res_eq;
    double *matrix_primal_part_1;
    double *matrix_primal;
    double *matrix_primal_part_2;
    double vector_primal_part;
    double vector_primal;  
    double dual;
    double primal_matrix_primal_part;
    double primal_matrix_primal;
    double complm_part; 
    double complm;   
    double sq_grad_x_part;
    double sq_nc_grad_x_part;
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
    double *solution_x;
    double *solution_lam;



    


    //Initiate Primal Communicator
    MPI_Comm_split(MPI_COMM_WORLD, !(rank < T1) ? MPI_UNDEFINED : 0, rank, &comm_1);
    if (rank < T1){
        MPI_Comm_size(comm_1, &size_1);
        MPI_Comm_rank(comm_1, &rank_1);

        //Initiate Primal Decision Variables:
        y = (double *)calloc(n_t, sizeof(double));
        x = (double *)calloc(n_t, sizeof(double));
        for (n_col = 0; n_col < n_t; n_col++){
            x[n_col] = 0.0;
        }
        if (!(rank_1)){
            v = (double *)calloc(1, sizeof(double));
            u = (double *)calloc(1, sizeof(double));
            u[0] = 0.0;
        }

        //Read in Problem Formulation:
        q_obj = (double *)calloc(n_t, sizeof(double));
        if (!(rank_1)){
        	c_obj = (double *)calloc(1, sizeof(double));
        	c_vec = (double **)calloc(M, sizeof(double *));
        	for (m = 0; m < M; m++){
        		c_vec[m] = (double *)calloc(1, sizeof(double));
        	}
            temp = (double *)calloc(N, sizeof(double));
        }
        P = (double ***)calloc(M, sizeof(double **));
        for (m = 0; m < M; m++){
            P[m] = (double **)calloc(N, sizeof(double *));
            for (n_row = 0; n_row < N; n_row++){
                P[m][n_row] = (double *)calloc(n_t, sizeof(double));
            }
        }       
        A = (double **)calloc(1, sizeof(double *));
        A[0] = (double *)calloc(n_t, sizeof(double));


        //Initiate Objective Function:
        //Linear Part:
        for (n_col = 0; n_col < n_t; n_col++){
        	q_obj[n_col] = -1.0;
        }      
        if (!(rank_1)){                                    
            c_obj[0] = (double)M;
        }
        //Initiate Linear Equality Constraints:
        if (!(rank_1)){
            sprintf(filename1, "%s/train_label", FOLDER);
            if(!(fp1 = fopen(filename1,"r"))){
                printf("Cannot open train_label file\n");
                exit(0);
            }
            for (n_col = 0; n_col < N; n_col++){
                fscanf(fp1, "%lf", &temp[n_col]);
            }
            fclose(fp1);
        }
        MPI_Scatter(temp, n_t, MPI_DOUBLE, A[0], n_t, MPI_DOUBLE, 0, comm_1);
        //Initiate Quadratic Inequality Constraints:
        for (m = 0; m < M; m++){           
            if (!(rank_1)){
                //Linear Part:
            	c_vec[m][0] = -1.0;
                //Quadratic Part:
                sprintf(filename1, "%s/%d/G_%d", FOLDER, M, m); //P_m is the file containing the quadratic part of constraint m
                if(!(fp1 = fopen(filename1,"r"))){
                    printf("Cannot open kernel file\n");
                    exit(0);
                }
                fscanf(fp1, "%d", &nzcnt);
            }
            for (n_row = 0; n_row < N; n_row++){
                if (!(rank_1)){
                    for (n_col = 0; n_col < N; n_col++){
                        fscanf(fp1, "%lf", &temp[n_col]);
                        temp[n_col] = temp[n_col]/NUM_DATA;
                    }
                }
                MPI_Scatter(temp, n_t, MPI_DOUBLE, P[m][n_row], n_t, MPI_DOUBLE, 0, comm_1);
            }
            if (!(rank_1)){
                fclose(fp1);
            }
        }       
        if (!(rank_1)){
printf("Initiation of Objective Function and Left-Hand Side of Constraints is OK!\n");
        }

        //Calculate Matrix Norm:        
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
        norm_A_part = 0.0;
        for (n_col = 0; n_col < n_t; n_col++){
        	norm_A_part += pow(A[0][n_col], (double)2);
        }
        MPI_Reduce(&norm_A_part, &norm_A, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        if (!(rank_1)){
            norm_PP = 0.0;
            for (m = 0; m < M; m++){
                norm_PP += norm_P[m];
                norm_P[m] = sqrt(norm_P[m]);
            }
            norm_PP = sqrt(norm_PP);

            norm_C = 0.0;
            for (m = 0; m < M; m++){
            	norm_C += pow(c_vec[m][0], (double)2);
            }
            norm_C = sqrt(norm_C);

            norm_A = sqrt(norm_A);
printf("Calculation of Matrix Norms is OK!\n");
        }

        //Initiate Our Algorithm Variables:
        grad_x = (double *)calloc(n_t, sizeof(double));
        nc_grad_x = (double *)calloc(n_t, sizeof(double));
        matrix_primal_part_1 = (double *)calloc(N, sizeof(double)); //(part of P)*(part of x)
        matrix_primal_part_2 = (double *)calloc(n_t, sizeof(double)); //part of P*x
        if (!(rank_1)){
        	grad_u = (double *)calloc(1, sizeof(double));
        	matrix_primal = (double *)calloc(N, sizeof(double));
            weight = (double *)calloc(6, sizeof(double));
            for (i = 1; i < 6; i++){
                weight[i] = 1.0;
            }        
            tol = (double *)calloc(2, sizeof(double));

            //Open File
            sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_%d_Tolerance", FOLDER, M, T1); //OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_T1_Tolerance is the file containing the tolerance of each iteration
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

        //Initiate Dual Decision Variables:
        mu = (double *)calloc(m_t, sizeof(double));
        lam = (double *)calloc(m_t, sizeof(double));
        for (m = 0; m < m_t; m++){
            lam[m] = 1.0;
        }
        if (!(rank_2)){
        	nu = (double *)calloc(1, sizeof(double));
        	gam = (double *)calloc(1, sizeof(double));
        	gam[0] = 1.0;
        }

        //Initiate Matrix Norm:
        norm_P_part = (double *)calloc(m_t, sizeof(double));
        if (!(rank_2)){
        	norm_P = (double *)calloc(M, sizeof(double));
        }

        //Initiate Our Algorithm Variables:
        res_ineq = (double *)calloc(m_t, sizeof(double));
        if (!(rank_2)){
        	res_eq = (double *)calloc(1, sizeof(double));
        }
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
        //Update Gradient and Residual:
        if (rank < T1){
        	//Initiate Gradient_x with q_obj:
        	for (n_col = 0; n_col < n_t; n_col++){
        		grad_x[n_col] =  q_obj[n_col];
        	}
        	//Calculate A*x:
        	vector_primal_part = 0.0;
        	for (n_col = 0; n_col < n_t; n_col++){
        		vector_primal_part += A[0][n_col]*x[n_col];
        	}
        	MPI_Reduce(&vector_primal_part, &vector_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	if (!(rank_1)){
        		//Initiate Gradient_u with c_obj:
        		grad_u[0] = c_obj[0];
        		//Send A*x:
        		MPI_Send(&vector_primal, 1, MPI_DOUBLE, T1, 233, MPI_COMM_WORLD);
        	}
        }
        else{
            complm_part = 0.0;
            if (!(rank_2)){
            	//Receive A*x and Update Equality Residual:
            	MPI_Recv(&res_eq[0], 1, MPI_DOUBLE, 0, 233, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for (m = 0; m < M; m++){
            //Send Dual:
        	if (rank == T1 + (m/m_t)){
        		MPI_Send(&lam[(m % m_t)], 1, MPI_DOUBLE, 0, 0*M + m, MPI_COMM_WORLD);
        	}
        	//Receive Dual:
        	if (!(rank)){
        		MPI_Recv(&dual, 1, MPI_DOUBLE, (T1 + m/m_t), 0*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        		//Update Gradient_u with lam_m*c_m:
        		grad_u[0] += dual*c_vec[m][0];
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
        		//Update Gradient_x with lam_m*(P_m*x):
        		for (n_col = 0; n_col < n_t; n_col++){
        			grad_x[n_col] += dual*matrix_primal_part_2[n_col];
        		}
        		//Calculate x*P_m*x:
        		primal_matrix_primal_part = 0.0;
        		for (n_col = 0; n_col < n_t; n_col++){
        			primal_matrix_primal_part += matrix_primal_part_2[n_col]*x[n_col];
        		}
        		MPI_Reduce(&primal_matrix_primal_part, &primal_matrix_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	}        	
        	if (!(rank)){
        		//Update with c_m*u:
        		primal_matrix_primal = 0.5*primal_matrix_primal + c_vec[m][0]*u[0];
        		//Send 0.5*x*P_m*x + c_m*u:
        		MPI_Send(&primal_matrix_primal, 1, MPI_DOUBLE, (T1 + m/m_t), 1*M + m, MPI_COMM_WORLD);
        	}
        	//Receive 0.5*x*P_m*x + c_m*u and Update Inequality Residual:
        	if (rank == T1 + (m/m_t)){
        		MPI_Recv(&res_ineq[(m % m_t)], 1, MPI_DOUBLE, 0, 1*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (res_ineq[(m % m_t)] >= 0){
                    complm_part += pow(lam[(m % m_t)]*res_ineq[(m % m_t)], (double)2);
                }
                else{
                    complm_part += pow(lam[(m % m_t)]*(-res_ineq[(m % m_t)]), (double)2);
                }
        	}
        }        
        //Send Dual:
        if (rank == T1){
        	MPI_Send(&gam[0], 1, MPI_DOUBLE, 0, 2333, MPI_COMM_WORLD);
        }
        //Receive Dual:
        if ((!rank)){
        	MPI_Recv(&dual, 1, MPI_DOUBLE, T1, 2333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < T1){
        	//Broadcast Dual:
        	MPI_Bcast(&dual, 1, MPI_DOUBLE, 0, comm_1);
        	sq_grad_x_part = 0.0;
        	sq_nc_grad_x_part = 0.0;
        	//Update Gradient_x with gam*A^T:
        	for (n_col = 0; n_col < n_t; n_col++){
        		grad_x[n_col] += dual*A[0][n_col];
        		sq_grad_x_part += pow(grad_x[n_col], (double)2);
        		//Calculate Normal Cone Gradient_x:
        		if (x[n_col] == 0){
        			if (grad_x[n_col] < 0.0){
        				nc_grad_x[n_col] = grad_x[n_col];
        			}
        			else{
        				nc_grad_x[n_col] = 0.0;
        			}
        		}
        		else{
                    if (x[n_col] == C){
                        if (grad_x[n_col] > 0.0){
                            nc_grad_x[n_col] = grad_x[n_col];
                        }
                        else{
                            nc_grad_x[n_col] = 0.0;
                        }
                    }
                    else{
                        nc_grad_x[n_col] = grad_x[n_col];
                    }
        		}
        		sq_nc_grad_x_part += pow(nc_grad_x[n_col], (double)2);
        	}
        	MPI_Reduce(&sq_grad_x_part, &a, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
            MPI_Reduce(&sq_nc_grad_x_part, &tol[0], 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        }        
        else{
            MPI_Reduce(&complm_part, &complm, 1, MPI_DOUBLE, MPI_SUM, 0, comm_2);
            if (rank == T1){
            	complm = complm + pow(res_eq[0], (double)2);
                MPI_Send(&complm, 1, MPI_DOUBLE, 0, 888, MPI_COMM_WORLD);
            }
        }
        if (!(rank)){
            MPI_Recv(&tol[1], 1, MPI_DOUBLE, T1, 888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            a = sqrt(a);
            
            //Calculate Tolerance:
            tol[0] = sqrt((tol[0] + pow(grad_u[0], (double)2))/(double)(N+1));
            tol[1] = sqrt(tol[1]/(double)(M+1));
            fprintf(fp1, "%.16lf\t%.16lf\n", tol[0], tol[1]);
//printf("%.16lf\t%.16lf\n", tol[0], tol[1]);
            if ((tol[0] < TOLERANCE_1) && (tol[1] < TOLERANCE_2)){
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
            MPI_Send(&epsilon[1], 1, MPI_DOUBLE, T1, 8888, MPI_COMM_WORLD);
        }

        //Update Stepsize:
        if (rank >= T1){
            //Rule 2:
            if (rank == T1){
                MPI_Recv(&epsilon[1], 1, MPI_DOUBLE, 0, 8888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Bcast(&epsilon[1], 1, MPI_DOUBLE, 0, comm_2);
            rho[1] = LARGE;
            for (m = 0; m < m_t; m++){
                if (res_ineq[m] < 0){
                    a = -res_ineq[m];
                }
                else{
                    a = res_ineq[m];
                }
                b = lam[m];
                if (norm_P_part[m] == 0.0){
                    c = epsilon[1]/(double)M;
                }
                else{
                    c = epsilon[1]/(M*norm_P_part[m]);
                }
                if (a > TOLERANCE_1){
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
                if (rho_2_tmp < rho[1]){
                    rho[1] = rho_2_tmp;
                }
            }
            MPI_Reduce(&rho[1], &rho_2_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm_2);
            if (rank == T1){
                MPI_Send(&rho_2_min, 1, MPI_DOUBLE, 0, 88888, MPI_COMM_WORLD);
            }
        }
        if (!(rank)){
        	//Rule 3:
        	b = 2*norm_x;
        	c = 2*epsilon[2]/norm_PP;
        	if (a > TOLERANCE_1){
        		rho[2] = (sqrt(pow(b, (double)2) + 4*a*c) - b)/(2*a);
        	}
        	else{
        		if (b == 0.0){
        			rho[2] = 2*epsilon[2];
        		}
        		else{
        			rho[2] = c/b;
        		}
        	}
            //Rule 5:
            if (norm_x == 0.0){
            	rho[3] = epsilon[3];
            }
            else{
            	rho[3] = epsilon[3]/(norm_x*norm_PP);
            }
            //Rule 6:
            if (norm_C == 0.0){
            	rho[4] = epsilon[4];
            }
            else{
            	rho[4] = epsilon[4]/norm_C;
            }
            //Rule 7:
            if (norm_A == 0.0){
            	rho[5] = epsilon[5];
            }
            else{
            	rho[5] = epsilon[5]/norm_A;
            }
            MPI_Recv(&rho[1], 1, MPI_DOUBLE, T1, 88888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
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
        	//x
        	for (n_col = 0; n_col < n_t; n_col++){
        		y[n_col] = x[n_col] - rho[0]*grad_x[n_col];
        		if (y[n_col] < 0.0){
        			y[n_col] = 0.0;
        		}
                if (y[n_col] > C){
                    y[n_col] = C;
                }
        	}
        	//u
        	if (!(rank_1)){
        		v[0] = u[0] - rho[0]*grad_u[0];
        	}
        }       
        //Dual
        else{
        	//Inequality
        	for (m = 0; m < m_t; m++){
        		mu[m] = lam[m] + rho[0]*res_ineq[m];
        		if (mu[m] < 0.0){
        			mu[m] = 0.0;
        		}
        	}
        	//Equality
        	if (!(rank_2)){
        		nu[0] = gam[0] + rho[0]*res_eq[0];
        	}
        }


        //Update Gradient and Residual:
        if (rank < T1){
        	//Initiate Gradient_x with q_obj:
        	for (n_col = 0; n_col < n_t; n_col++){
        		grad_x[n_col] =  q_obj[n_col];
        	}
        	//Calculate A*y:
        	vector_primal_part = 0.0;
        	for (n_col = 0; n_col < n_t; n_col++){
        		vector_primal_part += A[0][n_col]*y[n_col];
        	}
        	MPI_Reduce(&vector_primal_part, &vector_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	if (!(rank_1)){
        		//Initiate Gradient_u with c_obj:
        		grad_u[0] = c_obj[0];
        		//Send A*y:
        		MPI_Send(&vector_primal, 1, MPI_DOUBLE, T1, 23333, MPI_COMM_WORLD);
        	}
        }
        else{
            if (!(rank_2)){
            	//Receive A*y and Update Equality Residual:
            	MPI_Recv(&res_eq[0], 1, MPI_DOUBLE, 0, 23333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
        		//Update Gradient_u with mu_m*c_m:
        		grad_u[0] += dual*c_vec[m][0];
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
        		//Update Gradient_x with mu_m*(P_m*y):
        		for (n_col = 0; n_col < n_t; n_col++){
        			grad_x[n_col] += dual*matrix_primal_part_2[n_col];
        		}
        		//Calculate y*P_m*y:
        		primal_matrix_primal_part = 0.0;
        		for (n_col = 0; n_col < n_t; n_col++){
        			primal_matrix_primal_part += matrix_primal_part_2[n_col]*y[n_col];
        		}
        		MPI_Reduce(&primal_matrix_primal_part, &primal_matrix_primal, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
        	}        	
        	if (!(rank)){
        		//Update with c_m*v:
        		primal_matrix_primal = 0.5*primal_matrix_primal + c_vec[m][0]*v[0];
        		//Send 0.5*y*P_m*y + c_m*v:
        		MPI_Send(&primal_matrix_primal, 1, MPI_DOUBLE, (T1 + m/m_t), 3*M + m, MPI_COMM_WORLD);
        	}
        	//Receive 0.5*y*P_m*y + c_m*v and Update Inequality Residual:
        	if (rank == T1 + (m/m_t)){
        		MPI_Recv(&res_ineq[(m % m_t)], 1, MPI_DOUBLE, 0, 3*M + m, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	}
        }        
        //Send Dual:
        if (rank == T1){
        	MPI_Send(&nu[0], 1, MPI_DOUBLE, 0, 233333, MPI_COMM_WORLD);
        }
        //Receive Dual:
        if ((!rank)){
        	MPI_Recv(&dual, 1, MPI_DOUBLE, T1, 233333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < T1){
        	//Broadcast Dual:
        	MPI_Bcast(&dual, 1, MPI_DOUBLE, 0, comm_1);
        	//Update Gradient_x with nu*A^T:
        	for (n_col = 0; n_col < n_t; n_col++){
        		grad_x[n_col] += dual*A[0][n_col];
        	}
        }        


        //Update Primal and Dual Corrector:
        //Primal
        if (rank < T1){
        	//x
        	for (n_col = 0; n_col < n_t; n_col++){
        		x[n_col] += -rho[0]*grad_x[n_col];
        		if (x[n_col] < 0.0){
        			x[n_col] = 0.0;
        		}
                if (x[n_col] > C){
                    x[n_col] = C;
                }
        	}
        	//u
        	if (!(rank_1)){
        		u[0] += - rho[0]*grad_u[0];
        	}
        }       
        //Dual
        else{
        	//Inequality
        	for (m = 0; m < m_t; m++){
        		lam[m] += rho[0]*res_ineq[m];
        		if (lam[m] < 0.0){
        			lam[m] = 0.0;
        		}
        	}
        	//Equality
        	if (!(rank_2)){
        		gam[0] += rho[0]*res_eq[0];
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
    	objval_part = 0.0;
    	for (n_col = 0; n_col < n_t; n_col++){
    		objval_part += q_obj[n_col]*x[n_col];
    	}
    	MPI_Reduce(&objval_part, &objval, 1, MPI_DOUBLE, MPI_SUM, 0, comm_1);
    	if (!(rank_1)){
    		objval += c_obj[0]*u[0];
            solution_x = (double *)calloc(N, sizeof(double));
    	}
        MPI_Gather(x, n_t, MPI_DOUBLE, solution_x, n_t, MPI_DOUBLE, 0, comm_1);
    }
    else{
    	if (!(rank_2)){
    		solution_lam = (double *)calloc(M, sizeof(double));
    	}
    	MPI_Gather(lam, m_t, MPI_DOUBLE, solution_lam, m_t, MPI_DOUBLE, 0, comm_2);
    	if (!(rank_2)){
    		MPI_Send(solution_lam, M, MPI_DOUBLE, 0, 888888, MPI_COMM_WORLD);
            MPI_Send(gam, 1, MPI_DOUBLE, 0, 8888888, MPI_COMM_WORLD);
    	}
    }

    if (!(rank)){
    	solution_lam = (double *)calloc(M, sizeof(double));
    	MPI_Recv(solution_lam, M, MPI_DOUBLE, T1, 888888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        gam = (double *)calloc(1, sizeof(double));
        MPI_Recv(gam, 1, MPI_DOUBLE, T1, 8888888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	//Output Solution:
    	sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_%d", FOLDER, M, T1);
    	if(!(fp1 = fopen(filename1,"w"))){
    		printf("Cannot open this file\n");
    		exit(0);
    	}
    	fprintf(fp1, "Time use is %.16lf.\n", elapsed_time);
    	fprintf(fp1, "The Objective Function Value is %.16lf.\n", objval);
        fclose(fp1);

        sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_%d_Alpha", FOLDER, M, T1);
        if(!(fp1 = fopen(filename1,"w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (n_col = 0; n_col < N; n_col++){
            fprintf(fp1, "%.16lf\n", solution_x[n_col]);
        }
        fclose(fp1);

        sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_%d_U", FOLDER, M, T1);
        if(!(fp1 = fopen(filename1,"w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        fprintf(fp1, "%.16lf\n", u[0]);
        fclose(fp1);

        sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_%d_Lambda", FOLDER, M, T1);
        if(!(fp1 = fopen(filename1,"w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (m = 0; m < M; m++){
            fprintf(fp1, "%.16lf\n", solution_lam[m]);
        }
        fclose(fp1);

        sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_Our_Algorithm_Parallel_Core_%d_Gamma", FOLDER, M, T1);
        if(!(fp1 = fopen(filename1,"w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        fprintf(fp1, "%.16lf\n", gam[0]);
        fclose(fp1);
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


