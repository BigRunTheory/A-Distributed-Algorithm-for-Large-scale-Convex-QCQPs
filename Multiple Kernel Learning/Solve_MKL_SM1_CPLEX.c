#include <ilcplex/cplex.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


static char FOLDER[] = "2_norm";
//static char FOLDER[] = "HEPMASS";
//static char FOLDER[] = "Ionosphere";
//Model Parameters:
#define N 11009 //number of decision variables
#define M 3 //number of constraints
#define NUM_DATA 13760
#define NUM_TRAIN_DATA 11008
#define NUM_TEST_DATA 2752
#define TOLERANCE 0.001
#define C 2.0






int Solve_CPLEX (double *x, double* lam, double *objval, int *solstat, double *obj, int ccnt, int rcnt, int nzcnt, double *rhsval1, char *sense1, int *rmatbeg, int *rmatind, double *rmatval, int *linnzcnt, int *quadnzcnt, double *rhsval2, char sense2, int **linind, double **linval, int **quadrow, int **quadcol, double **quadval){
    int status = 0;
    CPXENVptr env = NULL;
    CPXLPptr lp  = NULL;
    char errmsg[1024];
    double *lb = (double *)calloc(N, sizeof(double));
    double *ub = (double *)calloc(N, sizeof(double));
    int n, m;
    for (n = 0; n < N-1; n++){
        lb[n] = 0.0;
        ub[n] = C;
    }
    lb[N-1] = -CPX_INFBOUND;
    ub[N-1] = CPX_INFBOUND;
 
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
        printf ("Failure to turn on preprocessing presolve.\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXsetintparam (env, CPXPARAM_Preprocessing_QCPDuals, CPX_QCPDUALS_FORCE);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failure to force the calculation of dual values.\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXsetdblparam (env, CPXPARAM_Barrier_QCPConvergeTol, TOLERANCE);
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

    status = CPXnewcols (env, lp, N, obj, lb, ub, NULL, NULL);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failed to populate problem.\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

    status = CPXaddrows (env, lp, ccnt, rcnt, nzcnt, rhsval1, sense1, rmatbeg, rmatind, rmatval, NULL, NULL);
//printf("%d\n", status);
    if ( status ) {
        printf ("Failed to add linear constraints.\n");
        CPXgeterrorstring (env, status, errmsg);
        printf ("%s\n", errmsg);
        goto TERMINATE;
    }

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




    //Calculate Lagrangian Multipliers for Quadratical Constraints
    int cols = CPXgetnumcols (env, lp);
    int numqs = CPXgetnumqconstrs (env, lp);
    double *dense = NULL;
    double *grad = NULL;
    int *slackind = NULL;
    double *slackval = NULL;
    int k;
    int surplus, len;
    int conetop;
    double zerotol = 1e-6;
    int ok;
    double maxabs;


    if ((dense = malloc(sizeof(*dense)*cols)) == NULL || (grad = malloc(sizeof(*grad)*cols)) == NULL || (slackind = malloc(sizeof(*slackind)*cols)) == NULL || (slackval = malloc (sizeof(*slackval)*cols)) == NULL){
        status = CPXERR_NO_MEMORY;
        goto TERMINATE;
    }

    for (m = 0; m < numqs; m++){
        /* Clear out dense slack and gradient vector. */
        for (n = 0; n < cols; n++){
            dense[n] = 0.0;
            grad[n] = 0.0;
        }

        /* Get dual slack vector and expand it to a dense vector. */
        status = CPXgetqconstrdslack (env, lp, m, &len, slackind, slackval, cols, &surplus);
        if ( status ) {
            printf ("Failed to get dual slack values for a quadratic constraint.\n");
            CPXgeterrorstring (env, status, errmsg);
            printf ("%s\n", errmsg);
            goto TERMINATE;
        }
        for (k = 0; k < len; k++){
            dense[slackind[k]] = slackval[k];
        }

        //Compute value of derivative at optimal solution.
        //The derivative of a quadratic constraint x^TQx + a^Tx + b <= 0 is Q^Tx + Qx + a.
        conetop = 1;
        for (k = 0; k < quadnzcnt[m]; k++){
            if (fabs(x[quadrow[m][k]]) > zerotol || fabs(x[quadcol[m][k]]) > zerotol){
                conetop = 0;
            }
            grad[quadcol[m][k]] += quadval[m][k]*x[quadrow[m][k]];
            grad[quadrow[m][k]] += quadval[m][k]*x[quadcol[m][k]];
        }

        for (k = 0; k < linnzcnt[m]; k++){
            grad[linind[m][k]] += linval[m][k];
            if (fabs(x[linind[m][k]]) > zerotol){
                conetop = 0;
            }
        }

        if ( conetop ) {
            fprintf (stderr, "#### WARNING: Cannot compute dual multipler at cone top!\n");
            status = CPXERR_BAD_ARGUMENT;
            goto TERMINATE;
        }
        else {
            ok = 0;
            maxabs = -1.0;

            //Compute lam[m] as slack/gradient.
            //We may have several indices to choose from and use the one with largest absolute value in the denominator.
            for (n = 0; n < cols; n++) {
                if (fabs(grad[n]) > zerotol) {
                    if (fabs(grad[n]) > maxabs ) {
                        lam[m] = fabs(dense[n]/grad[n]);
                        maxabs = fabs(grad[n]);
                    }
                    ok = 1;
                }
            }
            if ( !ok ) {
                //Dual slack is all 0. lam[m] can be anything, just set to 0.
                lam[m] = 0.0;
            }
        }
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
    int m, n_row, n_col, n_att;
    int z;
    double temp;
    FILE *fp1;
    char filename1[100];




    //CPLEX Parameters:
    double *x = (double *)calloc(N, sizeof(double)); //decision variables
    double *lam = (double *)calloc(M, sizeof(double)); //Lagrangian multipliers
    double objval; //objective function value
    int solstat;
    //Linear Objective
    double *obj = (double *)calloc(N, sizeof(double)); //linear objective function
    //Linear Equality Constraint
    int ccnt;
    int rcnt;
    int nzcnt;
    double *rhsval1;
    char *sense1;
    int *rmatbeg;
    int *rmatind;
    double *rmatval;
    //Quadratic Inequality Constraint
    int *linnzcnt = (int *)calloc(M, sizeof(int));
    int **linind = (int **)calloc(M, sizeof(int *));
    double **linval = (double **)calloc(M, sizeof(double *));
    int *quadnzcnt = (int *)calloc(M, sizeof(int));
    double *rhsval2 = (double *)calloc(M, sizeof(double));
    char sense2 = 'L';  
    int **quadrow = (int **)calloc(M, sizeof(int *));
    int **quadcol = (int **)calloc(M, sizeof(int *));
    double **quadval = (double **)calloc(M, sizeof(double *));
    
    
    //Initiate Objective Function:
    //Linear Part: 
    for (n_col = 0; n_col < N-1; n_col++){
        obj[n_col] = -1.0;
    }
    obj[N-1] = (double)(M);
printf("Initiation of Objective Function of QCQP is OK!\n");            
    
    //Initiate Constraints:
    //Linear Equality Constraint
    ccnt = 0;
    rcnt = 1;
    nzcnt = NUM_TRAIN_DATA;
    rhsval1 = (double *)calloc(rcnt, sizeof(double));  
    sense1 = (char *)calloc(rcnt, sizeof(char));
    rmatbeg = (int *)calloc(rcnt, sizeof(int));
    rhsval1[0] = 0.0;
    sense1[0] = 'E';
    rmatbeg[0] = 0;
    rmatind = (int *)calloc(nzcnt, sizeof(int));
    rmatval = (double *)calloc(nzcnt, sizeof(double));
    sprintf(filename1, "%s/train_label", FOLDER);
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
        rmatind[n_col] = n_col;
        fscanf(fp1, "%lf", &rmatval[n_col]);
    }
    fclose(fp1);   
    //Quadratic Inequality Constraint
    for (m = 0; m < M; m++){
        //Linear Part:
        linnzcnt[m] = 1;
        linind[m] = (int *)calloc(linnzcnt[m], sizeof(int));
        linval[m] = (double *)calloc(linnzcnt[m], sizeof(double));
        linind[m][0] = N-1;
        linval[m][0] = -1.0;
        //Quadratic Part:
        sprintf(filename1, "%s/%d/G_%d", FOLDER, M, m); //G_m is the file containing the quadratic part of constraint m
        if(!(fp1 = fopen(filename1,"r"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        fscanf(fp1, "%d", &quadnzcnt[m]);
        quadrow[m] = (int *)calloc(quadnzcnt[m], sizeof(int));
        quadcol[m] = (int *)calloc(quadnzcnt[m], sizeof(int));
        quadval[m] = (double *)calloc(quadnzcnt[m], sizeof(double));
        z = 0;
        for (n_col = 0; n_col < N-1; n_col++){
            for (n_row = 0; n_row < N-1; n_row++){
                fscanf(fp1, "%lf", &quadval[m][z]);
                quadval[m][z] = 0.5*quadval[m][z]/NUM_DATA;                
                if (quadval[m][z] != 0.0){
                    quadrow[m][z] = n_row;
                    quadcol[m][z] = n_col;
                    z++;
                }
            }
        }
        fclose(fp1);
        //Right-Hand Side:
        rhsval2[m] = 0.0;
    }
printf("Initiation of Constraints of QCQP is OK!\n");       
    
    
    //Solve QCQP using CPLEX:
    start = time(NULL);
    Solve_CPLEX (x, lam, &objval, &solstat, obj, ccnt, rcnt, nzcnt, rhsval1, sense1, rmatbeg, rmatind, rmatval, linnzcnt, quadnzcnt, rhsval2, sense2, linind, linval, quadrow, quadcol, quadval); 
    end = time(NULL);
    printf("Time use is %.16lf.\n", (double)(end-start));




    //Output the solutions
    sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_CPLEX", FOLDER, M);
    if(!(fp1 = fopen(filename1, "w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    fprintf(fp1, "Time use is %.16lf.\n", (double)(end-start));
    fprintf(fp1, "The Objective Function Value is %.16lf.\n", objval);
    fclose(fp1);

    sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_CPLEX_Alpha", FOLDER, M);
    if(!(fp1 = fopen(filename1,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
        fprintf(fp1, "%.16lf\n", x[n_col]);
    }
    fclose(fp1);

    sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM1_CPLEX_Lambda", FOLDER, M);
    if(!(fp1 = fopen(filename1, "w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (m = 0; m < M; m++){
        fprintf(fp1, "%.16lf\n", lam[m]);
    }
    fclose(fp1);
}
