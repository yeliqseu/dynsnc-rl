#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include "bats.h"

// batch sparse coding parameters
extern int    currbatch;
extern int    batchsent;              // number of packets sent from the current batch
extern int    s_count;    // number of recoded packets sent from the current batc
extern int    batchcount;   // number of received packets of the current receiving batch
extern int    dofcount;     // number of innovative packets contributed by the current receiving batch

// Reinforcement learning parameters and allocated memories for tables
extern double alpha;
extern double beta;
extern double gamma;
extern double epsilon;
extern int    maxbts;
extern double **Qtable;
extern double **visits;
extern int    nstate;
extern int    action[][2];
extern int    naction;
extern double pe;


int load_table(char *fname, double **table, int nrow, int ncol)
{
    int exist = 0;
    if( access( fname, F_OK ) != -1 ) {
        exist = 1;
        FILE *fp = fopen(fname, "r");
        for (int i=0; i<nrow; i++) {
            for (int j=0; j<ncol; j++) {
                fscanf(fp, "%lf \t", &(table[i][j]));
            }
            fscanf(fp, "\n");
        }
    } else {
        exist = 0;
    }
    return exist;
}

int save_table(char *fname, double **table, int nrow, int ncol)
{
    FILE *f = fopen(fname, "w");
    for (int i=0; i<nrow; i++) {
        for (int j=0; j<ncol; j++) {
            fprintf(f, "%.4f\t", table[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int derive_e_greedy_action(int r, double **table, int n_action)
{
    int act_id = -1;
    if (rand() % 100000 < epsilon * 100000) {
        // explore, pick a random action
        act_id = rand() % n_action;
    } else {
        // exploit
        // int state_id = r + c * (c + 1) / 2;
        int state_id = r;
        int tmp = rand() % n_action;
        int maxvalue = table[state_id][tmp];
        act_id = tmp;
        for (int i=0; i<n_action; i++) {
            if (table[state_id][i] > maxvalue) {
                maxvalue = table[state_id][i];
                act_id = i;
            }
            /*
            if (table[state_id][i] == maxvalue) {
                // break tie arbitrarily
                if (rand() % 2 == 0) {
                    act_id = i;
                }
            }
            */
        }
    }
    return act_id;
}

int derive_optimal_action(int r, double **table, int n_action)
{
    int act_id = -1;
    // exploit
    // int state_id = r + c * (c + 1) / 2;
    int state_id = r;
    int maxvalue = table[state_id][0];
    act_id = 0;
    for (int i=1; i<n_action; i++) {
        if (table[state_id][i] > maxvalue) {
            maxvalue = table[state_id][i];
            act_id = i;
        }
        if (table[state_id][i] == maxvalue) {
            // break tie arbitrarily
            if (rand() % 2 == 0) {
                act_id = i;
            }
        }
    }
    return act_id;
}