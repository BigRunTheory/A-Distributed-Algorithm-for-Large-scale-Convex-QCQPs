#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//Data Size
#define NUM_ATTRIBUTE 20
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


	//Input Data:
	double **data = (double **)calloc(NUM_DATA, sizeof(double *));
	for (n_row = 0; n_row < NUM_DATA; n_row++){
		data[n_row] = (double *)calloc((NUM_ATTRIBUTE + 1), sizeof(double));
	}
	//Positive Class:
	sprintf(filename1, "positive_class.csv");
	if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    n_row = 0;
    while (fgets(buf, 1024, fp1) != NULL){
    	if (!(sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &data[n_row][0], &data[n_row][1], &data[n_row][2], &data[n_row][3], &data[n_row][4], &data[n_row][5], &data[n_row][6], &data[n_row][7], &data[n_row][8], &data[n_row][9], &data[n_row][10], &data[n_row][11], &data[n_row][12], &data[n_row][13], &data[n_row][14], &data[n_row][15], &data[n_row][16], &data[n_row][17], &data[n_row][18], &data[n_row][19]))){
    		exit(0);
    	}
    	data[n_row][NUM_ATTRIBUTE] = 1.0;
    	n_row++;
    }
    fclose(fp1);
    //Negative Class:
    sprintf(filename1, "negative_class.csv");
	if(!(fp1 = fopen(filename1,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    n_row = NUM_DATA/2;
    while (fgets(buf, 1024, fp1) != NULL){
    	if (!(sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &data[n_row][0], &data[n_row][1], &data[n_row][2], &data[n_row][3], &data[n_row][4], &data[n_row][5], &data[n_row][6], &data[n_row][7], &data[n_row][8], &data[n_row][9], &data[n_row][10], &data[n_row][11], &data[n_row][12], &data[n_row][13], &data[n_row][14], &data[n_row][15], &data[n_row][16], &data[n_row][17], &data[n_row][18], &data[n_row][19]))){
    		exit(0);
    	}
    	data[n_row][NUM_ATTRIBUTE] = -1.0;
    	n_row++;
    }
    fclose(fp1);


    //Shuffle Data:
    int *index = (int *)calloc(NUM_DATA, sizeof(int));
    for (n_row = 0; n_row < NUM_DATA; n_row++){
    	index[n_row] = n_row;
    }
    int temp;
    srand (time(NULL));
    for (n_row = NUM_DATA - 1; n_row > 0; n_row--){
    	n_col = rand()%(n_row + 1);
    	temp = index[n_row];
    	index[n_row] = index[n_col];
    	index[n_col] = temp;
    }


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
            fprintf(fp1, "%.16lf\t", data[index[n_row]][n_att]);
        }
        fprintf(fp1, "\n");
        fprintf(fp2, "%.16lf\n", data[index[n_row]][NUM_ATTRIBUTE]);
    }
    fclose(fp1);
    fclose(fp2);


    //Output Test Data:
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
    for (n_row = NUM_TRAIN_DATA; n_row < NUM_DATA; n_row++){
    	for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
            fprintf(fp1, "%.16lf\t", data[index[n_row]][n_att]);
        }
        fprintf(fp1, "\n");
        fprintf(fp2, "%.16lf\n", data[index[n_row]][NUM_ATTRIBUTE]);
    }
    fclose(fp1);
    fclose(fp2);
}