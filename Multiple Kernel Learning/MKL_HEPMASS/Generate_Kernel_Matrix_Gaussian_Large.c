#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


//static char FOLDER[] = "2_norm";
static char FOLDER[] = "HEPMASS";
//static char FOLDER[] = "Ionosphere";
#define NUM_DATA 13760
#define NUM_ATTRIBUTE 28
#define NUM_TRAIN_DATA 11008
#define NUM_TEST_DATA 2752
#define M 9

/*Gaussian Kernel exp(-(x_1 - x_2)^T (x_1 - x_2)/(2*sigma^2)) with sigma = 1.0
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*pow(1.0, (double)2)));
                }
            }
*/

/*Linear Kernel x_1^T x_2
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += data[n_row][n_att]*data[n_col][n_att];
                    }
                }
            }
*/

/*Polynomial Kernel (x_1^T x_2 + c)^d with c = 1.0 and d = 2
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += data[n_row][n_att]*data[n_col][n_att];
                    }
                    K[n_row][n_col] = pow(K[n_row][n_col] + 1.0, (double)2);
                }
            }
*/

/*Exponential Kernel exp(beta*(x_1^T x_2)) with beta = 1.0
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += data[n_row][n_att]*data[n_col][n_att];
                    }
                    K[n_row][n_col] = exp(1.0*K[n_row][n_col]);
                }
            }
*/

/*Laplacian Kernel exp(-alpha*|x_1 - x_2|) with alpha = 1.0
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += fabs(data[n_row][n_att] - data[n_col][n_att]);
                    }
                    K[n_row][n_col] = exp(-(1.0/1.0)*K[n_row][n_col]);
                }
            }
*/






int main (int argc, char *argv[]){
    int n_row, n_col, n_att, m;
    int nz_cnt;
    FILE *fp;
    char filename[100];


    //Input Data:
    double **data = (double **)calloc(NUM_DATA, sizeof(double *));
    for (n_row = 0; n_row < NUM_DATA; n_row++){
        data[n_row] = (double *)calloc(NUM_ATTRIBUTE, sizeof(double));
    }

    sprintf(filename, "%s/train_attribute", FOLDER);
    if(!(fp = fopen(filename,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
        for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
            fscanf(fp, "%lf\t", &data[n_row][n_att]);
        }
        fscanf(fp, "\n");
    }
    fclose(fp);

    sprintf(filename, "%s/test_attribute", FOLDER);
    if(!(fp = fopen(filename,"r"))){
        printf("Cannot open this file\n");
        exit(0);
    }
    for (n_row = NUM_TRAIN_DATA; n_row < NUM_DATA; n_row++){
        for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
            fscanf(fp, "%lf\t", &data[n_row][n_att]);
        }
        fscanf(fp, "\n");
    }
    fclose(fp);

    
    //Calculate Kernel Matrices:
    double **K = (double **)calloc(NUM_DATA, sizeof(double *));
    double **NK = (double **)calloc(NUM_DATA, sizeof(double *));    
    for (n_row = 0; n_row < NUM_DATA; n_row++){
        K[n_row] = (double *)calloc(NUM_DATA, sizeof(double));
        NK[n_row] = (double *)calloc(NUM_DATA, sizeof(double));
    }
    for (m = 0; m < M; m++){
        if (m == 0){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*0.0001));
                }
            }
        }
        if (m == 1){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*0.001));
                }
            }
        }
        if (m == 2){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*0.01));
                }
            }
        }
        if (m == 3){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*0.1));
                }
            }
        }
        if (m == 4){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*1));
                }
            }
        }
        if (m == 5){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*10));
                }
            }
        }
        if (m == 6){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*100));
                }
            }
        }
        if (m == 7){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*1000));
                }
            }
        }
        if (m == 8){
            for (n_row = 0; n_row < NUM_DATA; n_row++){
                for (n_col = 0; n_col < NUM_DATA; n_col++){
                    K[n_row][n_col] = 0.0;
                    for (n_att = 0; n_att < NUM_ATTRIBUTE; n_att++){
                        K[n_row][n_col] += pow(data[n_row][n_att] - data[n_col][n_att], (double)2);
                    }
                    K[n_row][n_col] = exp(-K[n_row][n_col]/(2*10000));
                }
            }
        }

        //Normalize Kernel Matrix:       
        for (n_row = 0; n_row < NUM_DATA; n_row++){
            for (n_col = 0; n_col < NUM_DATA; n_col++){
                NK[n_row][n_col] = K[n_row][n_col]/(sqrt(K[n_row][n_row]*K[n_col][n_col]));
            }
        }

        //Calculate G(K,tr) = label[n_row]*label[n_col]*NK[n_row][n_col]
        double *train_label = (double *)calloc(NUM_TRAIN_DATA, sizeof(double));
        sprintf(filename, "%s/train_label", FOLDER);
        if(!(fp = fopen(filename,"r"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
            fscanf(fp, "%lf\n", &train_label[n_col]);
        }
        fclose(fp);

        double **G = (double **)calloc(NUM_TRAIN_DATA, sizeof(double *));
        for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
            G[n_row] = (double *)calloc(NUM_TRAIN_DATA, sizeof(double));
        }
        nz_cnt = 0;
        for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
            for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
                G[n_row][n_col] = train_label[n_row]*train_label[n_col]*NK[n_row][n_col];
                if (G[n_row][n_col] != 0.0){
                    nz_cnt++;
                }
            }
        }

        //Output G(K_tr):
        sprintf(filename, "%s/%d/G_%d", FOLDER, M, m);
        if(!(fp = fopen(filename, "w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        fprintf(fp, "%d\n", nz_cnt);
        for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
            for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
                fprintf(fp, "%.16lf\n", G[n_row][n_col]);
            }
        }
        fclose(fp);

        //Output NK_tr:
        sprintf(filename, "%s/%d/K_%d", FOLDER, M, m);
        if(!(fp = fopen(filename, "w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (n_row = 0; n_row < NUM_TRAIN_DATA; n_row++){
            for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
                fprintf(fp, "%.16lf\n", NK[n_row][n_col]);
            }
        }
        fclose(fp);

        //Output NK_tr_t:
        sprintf(filename, "%s/%d/T_%d", FOLDER, M, m);
        if(!(fp = fopen(filename, "w"))){
            printf("Cannot open this file\n");
            exit(0);
        }
        for (n_row = NUM_TRAIN_DATA; n_row < NUM_DATA; n_row++){
            for (n_col = 0; n_col < NUM_TRAIN_DATA; n_col++){
                fprintf(fp, "%.16lf\n", NK[n_row][n_col]);
            }
        }
        fclose(fp);
printf("Kernel Matrix %d is Done!\n", m+1);
    }
}
