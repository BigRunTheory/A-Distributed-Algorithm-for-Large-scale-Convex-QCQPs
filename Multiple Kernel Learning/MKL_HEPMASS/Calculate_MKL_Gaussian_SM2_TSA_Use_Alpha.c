#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


//static char FOLDER[] = "2_norm";
static char FOLDER[] = "HEPMASS";
//static char FOLDER[] = "Ionosphere";
//static char METHOD[] = "CPLEX";
static char METHOD[] = "Our_Algorithm_Parallel_Core_128";
//Model Parameters:
#define M 9 //number of constraints
#define NUM_DATA 13760
#define NUM_TRAIN_DATA 11008
#define NUM_TEST_DATA 2752
#define C 1.0






int main (int argc, char *argv[]){
    int m, n_row, n_col, n_att;
    int z;
    double temp;
    FILE *fp1;
    char filename1[100];


    //Read in Lambda:
    double *lam = (double *)calloc(M, sizeof(double));
    sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM2_%s_Lambda", FOLDER, M, METHOD);
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (m = 0; m < M; m++){
        fscanf(fp1, "%lf", &lam[m]);
    }
    fclose(fp1);
    //Read in Alpha:
    double *alpha = (double *)calloc(NUM_TRAIN_DATA, sizeof(double));
    sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM2_%s_Alpha", FOLDER, M, METHOD);
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
        fscanf(fp1, "%lf", &alpha[n_col]);
    }
    fclose(fp1);
    //Read in Train Label:
    double *train_label = (double *)calloc(NUM_TRAIN_DATA, sizeof(double));
    sprintf(filename1, "%s/train_label", FOLDER);
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
        fscanf(fp1, "%lf", &train_label[n_col]);
    }
    fclose(fp1); 


    //Read in Training-Block Kernel Matrix
    double **K = (double **)calloc(NUM_TRAIN_DATA, sizeof(double *));
    for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
        K[n_row] = (double *)calloc(NUM_TRAIN_DATA, sizeof(double));
    }
    for (m = 0; m < M; m++){
        sprintf(filename1, "%s/%d/K_%d", FOLDER, M, m);
        if(!(fp1 = fopen(filename1,"r"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
            for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
                fscanf(fp1, "%lf\n", &temp);
                K[n_row][n_col] += lam[m]*temp/NUM_DATA;
            }
        }
        fclose(fp1);
    }
    //Calculate Offset:
    double offset;
    int num_sv = 0;
    for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
        if (alpha[n_row] > 0.001){
            temp = 0.0;
            for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
                temp += K[n_row][n_col]*alpha[n_col]*train_label[n_col];
            }
            temp += alpha[n_row]*train_label[n_row]/C;
            temp = train_label[n_row] - temp;
printf("%.16lf\n", temp);
            offset += temp;
            num_sv++;
        }
    }
    offset = offset/num_sv;
printf("The offset is: %.16lf.\n", offset);  
    //Read in Mixed-Block Kernel Matrix
    double **T = (double **)calloc(NUM_TEST_DATA, sizeof(double *));
    for (n_row = 0; n_row < NUM_TEST_DATA; n_row++){
        T[n_row] = (double *)calloc(NUM_TRAIN_DATA, sizeof(double));
    }
    for (m = 0; m < M; m++){
        sprintf(filename1, "%s/%d/T_%d", FOLDER, M, m);
        if(!(fp1 = fopen(filename1,"r"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (n_row = 0; n_row < NUM_TEST_DATA; n_row++){
            for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
                fscanf(fp1, "%lf\n", &temp);
                T[n_row][n_col] += lam[m]*temp/NUM_DATA;
            }
        }
        fclose(fp1);
    }
    //Read in Test Label
    double *test_label = (double *)calloc(NUM_TEST_DATA, sizeof(double));
    sprintf(filename1, "%s/test_label", FOLDER);
    if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_col = 0; n_col < NUM_TEST_DATA; n_col++){
        fscanf(fp1, "%lf\n", &test_label[n_col]);
    }
    fclose(fp1);
    //Predict Test Label
    double *test_label_predict = (double *)calloc(NUM_TEST_DATA, sizeof(double));
    for (n_row = 0; n_row < NUM_TEST_DATA; n_row++){
        for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
            test_label_predict[n_row] += T[n_row][n_col]*alpha[n_col]*train_label[n_col];
        }
        test_label_predict[n_row] += offset;
        if (test_label_predict[n_row] < 0.0){
            test_label_predict[n_row] = -1.0;
        }
        else{
            test_label_predict[n_row] = 1.0;
        }
    }
    //Calculate TSA
    int count = 0;
    for (n_col = 0; n_col < NUM_TEST_DATA; n_col++){
        if (test_label_predict[n_col] == test_label[n_col]){
            count++;
        }
    }
    

    //Output the solutions
    sprintf(filename1, "%s/%d/OUT_Solve_MKL_SM2_%s_TSA", FOLDER, M, METHOD);
    if(!(fp1 = fopen(filename1, "w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    fprintf(fp1, "TSA is %.16lf.\n", (double)count/(double)NUM_TEST_DATA);
    fclose(fp1);
}

