#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Data Size
#define NUM_ATTRIBUTE 28
#define NUM_TRAIN_DATA 11008
#define NUM_TEST_DATA 2752
#define NUM_DATA 13760






int main (int argc, char *argv[]){
	int n_row, n_col, n_att;
	FILE *fp1;
	FILE *fp2;
	char filename1[100];
	char filename2[100];
	char buf[1024];


	//Input Train Data:
	double **train = (double **)calloc(NUM_TRAIN_DATA, sizeof(double *));
	for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
		train[n_row] = (double *)calloc((NUM_ATTRIBUTE + 1), sizeof(double));
	}
	sprintf(filename1, "train.csv");
	if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    n_row = 0;
    while (fgets(buf, 1024, fp1) != NULL){
    	if (!(sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &train[n_row][0], &train[n_row][1], &train[n_row][2], &train[n_row][3], &train[n_row][4], &train[n_row][5], &train[n_row][6], &train[n_row][7], &train[n_row][8], &train[n_row][9], &train[n_row][10], &train[n_row][11], &train[n_row][12], &train[n_row][13], &train[n_row][14], &train[n_row][15], &train[n_row][16], &train[n_row][17], &train[n_row][18], &train[n_row][19], &train[n_row][20], &train[n_row][21], &train[n_row][22], &train[n_row][23], &train[n_row][24], &train[n_row][25], &train[n_row][26], &train[n_row][27], &train[n_row][28]))){
    		exit(0);
    	}
    	n_row++;
    }
    fclose(fp1);
    //Output Train Data:
    sprintf(filename1, "train_attribute");
    if(!(fp1 = fopen(filename1,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    sprintf(filename2, "train_label");
    if(!(fp2 = fopen(filename2,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
    	for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
            fprintf(fp1, "%.16lf\t", train[n_row][n_att+1]);
        }
        fprintf(fp1, "\n");
        if (train[n_row][0] == 0.0){
        	fprintf(fp2, "%.16lf\n", -1.0);
        }
        else{
        	fprintf(fp2, "%.16lf\n", 1.0);
        }
    }
    fclose(fp1);
    fclose(fp2);



    //Input Test Data:
	double **test = (double **)calloc(NUM_TEST_DATA, sizeof(double *));
	for (n_row = 0; n_row < NUM_TEST_DATA; n_row++){
		test[n_row] = (double *)calloc((NUM_ATTRIBUTE + 1), sizeof(double));
	}
	sprintf(filename1, "test.csv");
	if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    n_row = 0;
    while (fgets(buf, 1024, fp1) != NULL){
    	if (!(sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &test[n_row][0], &test[n_row][1], &test[n_row][2], &test[n_row][3], &test[n_row][4], &test[n_row][5], &test[n_row][6], &test[n_row][7], &test[n_row][8], &test[n_row][9], &test[n_row][10], &test[n_row][11], &test[n_row][12], &test[n_row][13], &test[n_row][14], &test[n_row][15], &test[n_row][16], &test[n_row][17], &test[n_row][18], &test[n_row][19], &test[n_row][20], &test[n_row][21], &test[n_row][22], &test[n_row][23], &test[n_row][24], &test[n_row][25], &test[n_row][26], &test[n_row][27], &test[n_row][28]))){
    		exit(0);
    	}
    	n_row++;
    }
    fclose(fp1);
    //Output Train Data:
    sprintf(filename1, "test_attribute");
    if(!(fp1 = fopen(filename1,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    sprintf(filename2, "test_label");
    if(!(fp2 = fopen(filename2,"w"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_row = 0; n_row < NUM_TEST_DATA; n_row++){
    	for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
            fprintf(fp1, "%.16lf\t", test[n_row][n_att+1]);
        }
        fprintf(fp1, "\n");
        if (test[n_row][0] == 0.0){
        	fprintf(fp2, "%.16lf\n", -1.0);
        }
        else{
        	fprintf(fp2, "%.16lf\n", 1.0);
        }
    }
    fclose(fp1);
    fclose(fp2);
}