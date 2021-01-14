#include <ilcplex/cplex.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Model Parameters:
#define N 16384 //number of decision variables
#define M 1 //number of constraints
#define CONDITION_NUMBER 1.25






int Solve_CPLEX (double *x, double *objval, int *solstat, double *obj, int qnzcnt, int *qmatbeg, int *qmatcnt, int *qmatind, double *qmatval, int *linnzcnt, int *quadnzcnt, double *rhsval2, char sense2, int **linind, double **linval, int **quadrow, int **quadcol, double **quadval){
    int status = 0;
	CPXENVptr env = NULL;
    CPXLPptr lp  = NULL;
    char errmsg[1024];
    double *lb = (double *)calloc(N, sizeof(double));
    int n;
    for (n = 0; n < N; n++){
        lb[n] = -CPX_INFBOUND;
    }
 
    env = CPXopenCPLEX (&status);
//printf("%d\n", status);
    if ( env == NULL ) {
       printf ("Could not open CPLEX environment.\n");
       CPXgeterrorstring (env, status, errmsg);
       printf ("%s\n", errmsg);
       goto TERMINATE;
    }

    status = CPXsetintparam (env, CPXPARAM_Threads, 20);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failure to set maximal number of parallel threads.\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
//printf("%d\n", status);
    if ( status ) {
    	printf ("Failure to turn on screen indicator.\n");
    	CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
    	goto TERMINATE;
    }

    status = CPXsetintparam (env, CPXPARAM_Preprocessing_Presolve, CPX_ON);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failure to turn off preprocessing presolve.\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXsetdblparam (env, CPXPARAM_Barrier_QCPConvergeTol, 0.001);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failure to set the tolerance on complementarity for convergence in quadratically constrained problems (QCPs).\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    lp = CPXcreateprob (env, &status, "QCQP");
//printf("%d\n", status);
    if ( lp == NULL ) {
        printf ("Failed to create LP.\n");
        CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    CPXchgobjsen (env, lp, CPX_MIN);

    status = CPXnewcols (env, lp, N, obj, lb, NULL, NULL, NULL);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failed to populate problem.\n");
        CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failed to add quadratic objective function.\n");
        CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    /*
    status = CPXaddrows (env, lp, ccnt, rcnt, nzcnt, rhsval1, sense1, rmatbeg, rmatind, rmatval, NULL, NULL);
printf("%d\n", status);
    if ( status ) {
        printf ("Failed to add linear constraints.\n");
        CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
        goto TERMINATE;
    }
    */

    int m;
    for (m = 0; m < M; m++){
    	status = CPXaddqconstr (env, lp, linnzcnt[m], quadnzcnt[m], rhsval2[m], sense2, linind[m], linval[m], quadrow[m], quadcol[m], quadval[m], NULL);
//printf("%d\n", status);
        if ( status ) {
        	printf ("Failed to add quadratic constraints.\n");
        	CPXgeterrorstring (env, status, errmsg);
        	printf ("%s\n", errmsg);
        	goto TERMINATE;
        }
    }

    status = CPXbaropt (env, lp);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failed to optimize QCQP.\n");
        CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXsolution (env, lp, solstat, objval, x, NULL, NULL, NULL);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failed to obtain solution.\n");
        CPXgeterrorstring (env, status, errmsg);
    	printf ("%s\n", errmsg);
        goto TERMINATE;
    }

TERMINATE:

    if ( lp != NULL ) {
        status = CPXfreeprob (env, &lp);
        if ( status ) {
            printf ("CPXfreeprob failed.\n");
            CPXgeterrorstring (env, status, errmsg);
            printf ("%s\n", errmsg);
        }
    }

    if ( env != NULL ) {
        status = CPXcloseCPLEX (&env);
        if ( status ) {
            char errmsg[1024];
            printf ("Could not close CPLEX environment.\n");
            CPXgeterrorstring (env, status, errmsg);
            printf ("%s\n", errmsg);
        }
    }
	
	return (status);
}






int main (int argc, char *argv[]){
	time_t start, end;
    int m, n_row, n_col;
    int z;


    //CPLEX Parameters:
    double *x = (double *)calloc(N, sizeof(double)); //decision variables
    double objval; //objective function value
    int solstat;
    double *obj = (double *)calloc(N, sizeof(double)); //linear objective function
    int qnzcnt; //total number of non-zero entries in quadratic objective function
    int *qmatbeg = (int *)calloc(N, sizeof(int)); //for each column j, the beginning index of the first non-zero entry
    int *qmatcnt = (int *)calloc(N, sizeof(int)); //number of non-zero entries in each column j in quadratic objective function
    int *qmatind; //the row number of the corresponding non-zero entry
    double *qmatval; //the value of the corresponding non-zero entry
    /*
    int ccnt;
    int rcnt;
    int nzcnt;
    double *rhsval1;
    char *sense1;
    int *rmatbeg;
    int *rmatind;
    double *rmatval;
    */
    int *linnzcnt = (int *)calloc(M, sizeof(int));
    int *quadnzcnt = (int *)calloc(M, sizeof(int));
    double *rhsval2 = (double *)calloc(M, sizeof(double));
    char sense2 = 'L';
    int **linind = (int **)calloc(M, sizeof(int *));
    double **linval = (double **)calloc(M, sizeof(double *));
    int **quadrow = (int **)calloc(M, sizeof(int *));
    int **quadcol = (int **)calloc(M, sizeof(int *));
    double **quadval = (double **)calloc(M, sizeof(double *));


    FILE *fp1;
    char filename1[100];
	
	
	//Initiate Objective Function:
	//Linear Part: 
    sprintf(filename1, "%d_%d/Condition_Number_%.2f/q_obj_%d_%d", M, N, CONDITION_NUMBER, M, N); //q_obj_M_N is the file containing the linear objective function
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < N; n_col++){
    	fscanf(fp1, "%lf", &obj[n_col]);
//printf("%.16lf\n", obj[n_col]);
    }
    fclose(fp1);

    //Quadratic Part:
    sprintf(filename1, "%d_%d/Condition_Number_%.2f/P_obj_%d_%d", M, N, CONDITION_NUMBER, M, N); //P_obj_M_N is the file containing the quadratic objective function
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    fscanf(fp1, "%d", &qnzcnt);
    qmatind = (int *)calloc(qnzcnt, sizeof(int));
    qmatval = (double *)calloc(qnzcnt, sizeof(double));
    z = 0;
    for (n_col = 0; n_col < N; n_col++){
    	qmatcnt[n_col] = 0;
    	qmatbeg[n_col] = z;
    	for (n_row = 0; n_row < N; n_row++){
    		fscanf(fp1, "%lf", &qmatval[z]);
//printf("%.16lf\t", qmatval[z]);
    		if (qmatval[z] != 0.0){
    			qmatind[z] = n_row;
//printf("%d\t", qmatind[z]);
                qmatcnt[n_col]++;
    			z++;
    		}
    	}
//printf("\n");
    }
    fclose(fp1);   
printf("Initiation of Objective Function is OK!\n");		
	
	
	//Initiate Constraints:
    for (m = 0; m < M; m++){
    	//Linear Part:
    	sprintf(filename1, "%d_%d/Condition_Number_%.2f/q_%d_%d_%d", M, N, CONDITION_NUMBER, m, M, N); //q_m_M_N is the file containing the linear part of constraint m
    	if(!(fp1 = fopen(filename1,"r"))){
    		printf("Cannot open this file\n");
    		exit(0);
    	}
    	fscanf(fp1, "%d", &linnzcnt[m]);
    	linind[m] = (int *)calloc(linnzcnt[m], sizeof(int));
    	linval[m] = (double *)calloc(linnzcnt[m], sizeof(double));
    	z = 0;
    	for (n_col = 0; n_col < N; n_col++){
    		fscanf(fp1, "%lf", &linval[m][z]);
//printf("%.16lf\t", linval[m][z]);
    		if (linval[m][z] != 0.0){
    			linind[m][z] = n_col;
//printf("%d\n", linind[m][z]);
    			z++;
    		}
    	}
    	fclose(fp1);

    	//Quadratic Part:
    	sprintf(filename1, "%d_%d/Condition_Number_%.2f/P_%d_%d_%d", M, N, CONDITION_NUMBER, m, M, N); //P_m_M_N is the file containing the quadratic part of constraint m
    	if(!(fp1 = fopen(filename1,"r"))){
    		printf("Cannot open this file\n");
    		exit(0);
    	}
    	fscanf(fp1, "%d", &quadnzcnt[m]);
    	quadrow[m] = (int *)calloc(quadnzcnt[m], sizeof(int));
    	quadcol[m] = (int *)calloc(quadnzcnt[m], sizeof(int));
    	quadval[m] = (double *)calloc(quadnzcnt[m], sizeof(double));
    	z = 0;
    	for (n_col = 0; n_col < N; n_col++){
    		for (n_row = 0; n_row < N; n_row++){
    			fscanf(fp1, "%lf", &quadval[m][z]);
                quadval[m][z] = 0.5*quadval[m][z];
//printf("%.16lf\t", quadval[m][z]);                
    			if (quadval[m][z] != 0.0){
    				quadrow[m][z] = n_row;
//printf("%d ", quadrow[m][z]);
    				quadcol[m][z] = n_col;
//printf("%d\t", quadcol[m][z]);
    				z++;
    			}
    		}
//printf("\n");
    	}
    	fclose(fp1);
    }

    //Right-Hand Side:
    sprintf(filename1, "%d_%d/Condition_Number_%.2f/r_%d_%d", M, N, CONDITION_NUMBER, M, N); //r_M_N is the file containing all the scalars of m constraints
    if(!(fp1 = fopen(filename1,"r"))){
    	printf("Cannot open this file\n");
    	exit(0);
    }
    for (m = 0; m < M; m++){
    	fscanf(fp1, "%lf", &rhsval2[m]);
    	rhsval2[m] = -rhsval2[m];
    }
    fclose(fp1);
printf("Initiation of Constraints is OK!\n");		
	
	
	//Solve the problem using CPLEX:
	start = time(NULL);
	Solve_CPLEX (x, &objval, &solstat, obj, qnzcnt, qmatbeg, qmatcnt, qmatind, qmatval, linnzcnt, quadnzcnt, rhsval2, sense2, linind, linval, quadrow, quadcol, quadval);	
	end = time(NULL);
    printf("Time use is %.16lf.\n", (double)(end-start));	
	
	
	//Output the solutions
	sprintf(filename1, "%d_%d/Condition_Number_%.2f/OUT_Solve_QCQP_CPLEX_%d_%d", M, N, CONDITION_NUMBER, M, N);
	if(!(fp1 = fopen(filename1, "w"))){
		printf("Cannot open this file\n");
		exit(0);
	}
	fprintf(fp1, "Time use is %.16lf.\n", (double)(end-start));
	fprintf(fp1, "The objective function value is %.16lf.\n", objval);
	fprintf(fp1, "The solution is:\n");
	for (n_col = 0; n_col < N; n_col++){
		fprintf(fp1, "%.16lf\n", x[n_col]);
	}
	fclose(fp1);
}
