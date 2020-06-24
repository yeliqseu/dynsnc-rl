/*
 * A BATS decoder which succeeds as soon as a sufficient number of encoded packets are received. The decoder can determine the innovativeness of coded packets upon reception. To reduce the decoding cost of the straightforward on-the-fly Gaussian elimination, the following decoder pays extra attention to the vector sparseness during  processing each received packet.  
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "bats.h"
#include "galois.h"

static int process_vector_inbatch(struct bats_decoder_ref *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message);
static int process_vector(struct bats_decoder_ref *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message);
static int apply_parity_check_matrix(struct bats_decoder_ref *dec_ctx);
static void back_substitution(struct bats_decoder_ref *dec_ctx);
static void bats_free_decoder_currbatch(struct bats_decoder_ref *dec_ctx);

extern void init_genrand(unsigned long s);

static int applying_precode = 0;

extern int batchcount;   // number of received packets of the current receiving batch
extern int dofcount;     // number of innovative packets contributed by the current receiving batch
extern int currbatch;

struct bats_decoder_ref *bats_create_decoder_ref(BATSparam *param)
{
    static char fname[] = "bats_create_decoder_ref";
    struct bats_decoder_ref *dctx = malloc(sizeof(struct bats_decoder_ref));
    dctx->param = param;
    // replicate the precode bipartite graph at the decoder side
    // init mt19937 RNG
    init_genrand(param->seed);
    // create bipartite graph
    if (param->cnum != 0) {
        if ( (dctx->graph = malloc(sizeof(BP_graph))) == NULL ) {
            fprintf(stderr, "%s: malloc BP_graph\n", fname);
            return NULL;
        }
        if (create_bipartite_graph(dctx->graph, param->snum, param->cnum) < 0)
            return NULL;
    }

    // setup other decoder state information
    dctx->received = 0;
    dctx->overhead = 0;
    dctx->DoF = 0;
    dctx->seen = calloc(param->snum+param->cnum, sizeof(int));
    dctx->covered = 0;
    if (param->cnum != 0) {
        for (int i=0; i<param->snum+param->cnum; i++) {
            dctx->seen[i] = 0;
        }
        // dctx->covered = param->snum+param->cnum;
        dctx->covered = 0;
    }
    dctx->de_precode = 0;
    dctx->finished = 0;
    dctx->operations = 0;

    dctx->row = (struct row_vector **) calloc(param->snum+param->cnum, sizeof(struct row_vector *));
    if (dctx->row == NULL) {
        fprintf(stderr, "%s: calloc dctx->row failed\n", fname);
        goto AllocError;
    }

    dctx->message = calloc(param->snum+param->cnum, sizeof(GF_ELEMENT*));
    if (dctx->message == NULL) {
        fprintf(stderr, "%s: calloc dctx->message failed\n", fname);
        goto AllocError;
    }
    for (int i=0; i<param->snum+param->cnum; i++) {
        dctx->message[i] = calloc(param->pktsize, sizeof(GF_ELEMENT));
        if (dctx->message[i] == NULL) {
            fprintf(stderr, "%s: calloc dctx->message[%d] failed\n", fname, i);
            goto AllocError;
        }
    }
    dctx->pp = calloc(param->snum+param->cnum, sizeof(GF_ELEMENT*));

    // small temporaray matrix
    dctx->currbid = -1;
    dctx->batch_row = NULL;
    dctx->batch_msg = NULL;
    // memory of batch_row and batch_msg are allocated when the batch is first seen.
    return dctx;

AllocError:
    bats_free_decoder_ref(dctx);
    return NULL;
}

void bats_free_decoder_ref(struct bats_decoder_ref *decoder)
{
    // Yes, I know there is memory leak here, but I don't care!
    int i, j;
    int snum = decoder->param->snum;
    int cnum = decoder->param->cnum;
    for (i=0; i<snum+cnum; i++) {
        if (decoder->row[i] != NULL) {
            free(decoder->row[i]->elem);
            free(decoder->row[i]);
            decoder->row[i] = NULL;
        }
    }
    free(decoder->row);
    for (i=0; i<snum+cnum; i++) {
        if (decoder->message[i] != NULL) {
            free(decoder->message[i]);
            decoder->message[i] = NULL;
        }
    }
    free(decoder->message);
    for (i=0; i<snum+cnum; i++) {
        if (decoder->pp[i] != NULL) {
            free(decoder->pp[i]);
            decoder->pp[i] = NULL;
        }
    }
    free(decoder->pp);
    if (decoder->batch_row != NULL) {
        bats_free_decoder_currbatch(decoder);
    }
    free(decoder);
    return;
}

static void bats_free_decoder_currbatch(struct bats_decoder_ref *dec_ctx)
{
    int i;
    for (i=0; i<dec_ctx->currpnum; i++) {
        if (dec_ctx->batch_row[i] != NULL) {
            free(dec_ctx->batch_row[i]->elem);
            free(dec_ctx->batch_row[i]);
            dec_ctx->batch_row[i] = NULL;
        }
        if (dec_ctx->batch_msg[i] != NULL) {
            free(dec_ctx->batch_msg[i]);
            dec_ctx->batch_msg[i] = NULL;
        }
    }
    free(dec_ctx->batch_row);
    dec_ctx->batch_row = NULL;
    free(dec_ctx->batch_msg);
    dec_ctx->batch_msg = NULL;
}

// process received BATS packet
int bats_process_packet_ref(struct bats_decoder_ref *dec_ctx, BATSpacket *pkt)
{
    static char fname[] = "process_packet_ref";
    if (pkt == NULL)
        return 0;

    int curr_DoF = dec_ctx->DoF;

    if (pkt->batchid != dec_ctx->currbid) {
        // need to clean up memory of the last batch
        // but I'm not going to do this for now. I don't care 
        // memory leak at this moment.
        if (dec_ctx->batch_row != NULL) {
            bats_free_decoder_currbatch(dec_ctx);
        }
        dec_ctx->currbid = pkt->batchid;
        dec_ctx->currpnum = pkt->degree;
        dec_ctx->batch_row = (struct row_vector **) calloc(dec_ctx->currpnum, sizeof(struct row_vector *));
        dec_ctx->batch_msg = calloc(dec_ctx->currpnum, sizeof(GF_ELEMENT*));
    }

    batchcount += 1;
    // update the seen list
    for (int i=0; i<pkt->degree; i++) {
        if (dec_ctx->seen[pkt->pktid[i]] == 0) {
            dec_ctx->seen[pkt->pktid[i]] = 1;
            dec_ctx->covered += 1;
        }
    }

    dec_ctx->received += 1;
    dec_ctx->overhead += 1;

    // Process the received packet in the batch first
    /*
    int pivot_ib = process_vector_inbatch(dec_ctx, pkt->coes, pkt->syms);

    if (pivot_ib == -1) {
        // not an innovative BEV, return
        return;
    }
    */

    // process the packet against the global decoding matrix
    int i, j, k;
    GF_ELEMENT quotient;

    int pktsize = dec_ctx->param->pktsize;
    int numpp = dec_ctx->param->snum + dec_ctx->param->cnum;


    // Transform BATS packet's encoding vector to full length
    GF_ELEMENT *ces = calloc(numpp, sizeof(GF_ELEMENT));
    if (ces == NULL)
        fprintf(stderr, "%s: calloc ces failed\n", fname);
    for (i=0; i<pkt->degree; i++)
        ces[pkt->pktid[i]] = pkt->coes[i];

    // Process full-length encoding vector against decoding matrix
    int lastDoF = dec_ctx->DoF;
    int pivot = process_vector(dec_ctx, ces, pkt->syms);

    int newDoF = pivot >=0 ? 1 : 0;
    printf("[Batch %d] Received-DoF: %d New-DoF: %d\n", currbatch, curr_DoF, newDoF);

    free(ces);
    ces = NULL;
    // Apply parity-check vectors
    if (dec_ctx->DoF == dec_ctx->param->snum && dec_ctx->param->cnum != 0) {
        dec_ctx->de_precode = 1;    /*Mark de_precode before applying precode matrix*/
        
        applying_precode = 1;
        int missing_DoF = apply_parity_check_matrix(dec_ctx);
        applying_precode = 0;
        printf("After applying the parity-check matrix, %d DoF are missing.\n", missing_DoF);
        dec_ctx->DoF = numpp - missing_DoF;
        // double check how many packets are not seen yet after applying precode
        int total_covered = 0;
        for (int i=0; i<numpp; i++) {
            if (dec_ctx->row[i] != NULL) {
                total_covered += 1;
                dec_ctx->seen[i] = 1;       // mark seen list as well
                continue;                   // covered for sure
            }
            else {
                // need to check if above rows have the packet covered. Example:
                // 1 x x 0
                // 0 0 x 0
                // 0 0 x 0
                // 0 0 0 0
                // The second packet IS covered and the last one IS NOT.
                int covered = 0;
                for (int j=0; j<i; j++) {
                    if (dec_ctx->row[j] == NULL || dec_ctx->row[j]->len < i-j+1)
                        continue;
                    else {
                        if (dec_ctx->row[j]->elem[i-j] != 0) {
                            covered = 1;
                            dec_ctx->seen[i] = 1;
                            break;          // packet i is covered in the j-th row
                        }
                    }
                }
                if (covered)
                    total_covered += 1;
            }
        }
        dec_ctx->covered = total_covered;
        printf("After applying precode, %d packets are covered\n", total_covered);
    }

    if (dec_ctx->DoF == dec_ctx->param->snum + dec_ctx->param->cnum) {
        back_substitution(dec_ctx);
        printf("Received %d packets from batch %d and %d are innovative\n", batchcount, currbatch, dofcount);
    }
}

static int process_vector_inbatch(struct bats_decoder_ref *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message)
{
    static char fname[] = "process_vector";
    int i, j, k;
    int pivot = -1;
    int pivotfound = 0;
    GF_ELEMENT quotient;

    int pktsize = dec_ctx->param->pktsize;
    int numpp   = dec_ctx->currpnum;

    int rowop = 0;
    for (i=0; i<numpp; i++) {
        if (vector[i] != 0) {
            if (dec_ctx->batch_row[i] != NULL) {
                /* There is a valid row saved for pivot-i, process against it */
                quotient = galois_divide(vector[i], dec_ctx->batch_row[i]->elem[0]);
                galois_multiply_add_region(&(vector[i]), dec_ctx->batch_row[i]->elem, quotient, dec_ctx->batch_row[i]->len);
                galois_multiply_add_region(message, dec_ctx->batch_msg[i], quotient, pktsize);
                dec_ctx->operations += 1 + dec_ctx->batch_row[i]->len + pktsize;
                rowop += 1;
            } else {
                pivotfound = 1;
                pivot = i;
                break;
            }
        }
    }

    if (pivotfound == 1) {
        /* Save it to the corresponding row */
        dec_ctx->batch_row[pivot] = (struct row_vector*) malloc(sizeof(struct row_vector));
        if (dec_ctx->batch_row[pivot] == NULL)
            fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, pivot);
        int len = numpp - pivot;
        dec_ctx->batch_row[pivot]->len = len;

        dec_ctx->batch_row[pivot]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
        if (dec_ctx->batch_row[pivot]->elem == NULL)
            fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, pivot);
        memcpy(dec_ctx->batch_row[pivot]->elem, &(vector[pivot]), len*sizeof(GF_ELEMENT));

        dec_ctx->batch_msg[pivot] = (GF_ELEMENT *) calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dec_ctx->batch_msg[pivot], message,  pktsize*sizeof(GF_ELEMENT));
    }
    return pivot;
}

static int process_vector(struct bats_decoder_ref *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message)
{
    static char fname[] = "process_vector";
    int i, j, k;
    int pivot = -1;
    int pivotfound = 0;
    GF_ELEMENT quotient;

    int pktsize = dec_ctx->param->pktsize;
    int numpp   = dec_ctx->param->snum + dec_ctx->param->cnum;

    int rowop = 0;
    for (i=0; i<numpp; i++) {
        if (vector[i] != 0) {
            // calculate the density of the vector
            int density_c = 0;
            for (j=i; j<numpp; j++) {
                if (vector[j] != 0)
                    density_c += 1;
            }

            // calculate the length between the leading element and the last nonzero
            int nzlen_c = 0;  
            for (j=numpp-1; j>=0; j--) {
                if (vector[j] != 0) {
                    nzlen_c = j - i + 1;
                    break;
                }
            }

            if (dec_ctx->row[i] != NULL) {
                // There is a valid row saved for pivot-i, process against it
                // But swap if the vector is sparser than the stored one before processing.
                int density = 0;
                for (j=0; j<dec_ctx->row[i]->len; j++) {
                    if (dec_ctx->row[i]->elem[j] != 0)
                        density += 1;
                }

                int nzlen = 0;
                for (j=dec_ctx->row[i]->len-1; j>=0; j--) {
                    if (dec_ctx->row[i]->elem[j] != 0) {
                        nzlen = j - i + 1;
                        break;
                    }
                }

                // perform swap
                if (density_c < density) {
                    for (j=0; j<dec_ctx->row[i]->len; j++) {
                        GF_ELEMENT temp = dec_ctx->row[i]->elem[j];
                        dec_ctx->row[i]->elem[j] = vector[i+j];
                        vector[i+j] = temp;
                    }
                    //GF_ELEMENT *temp = dec_ctx->message[i];
                    //dec_ctx->message[i] = message;
                    // message = temp;
                    for (j=0; j<pktsize; j++) {
                        GF_ELEMENT temp = dec_ctx->message[i][j];
                        dec_ctx->message[i][j] = message[j];
                        message[j] = temp;
                    }
                }

                quotient = galois_divide(vector[i], dec_ctx->row[i]->elem[0]);
                galois_multiply_add_region(&(vector[i]), dec_ctx->row[i]->elem, quotient, dec_ctx->row[i]->len);
                galois_multiply_add_region(message, dec_ctx->message[i], quotient, pktsize);
                dec_ctx->operations += 1 + dec_ctx->row[i]->len + pktsize;
                rowop += 1;
            } else {
                pivotfound = 1;
                pivot = i;
                break;
            }
        }
    }

    if (pivotfound == 1) {
        /* Save it to the corresponding row */
        dec_ctx->row[pivot] = (struct row_vector*) malloc(sizeof(struct row_vector));
        if (dec_ctx->row[pivot] == NULL)
            fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, pivot);
        int len = numpp - pivot;
        dec_ctx->row[pivot]->len = len;
        dec_ctx->row[pivot]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
        if (dec_ctx->row[pivot]->elem == NULL)
            fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, pivot);
        memcpy(dec_ctx->row[pivot]->elem, &(vector[pivot]), len*sizeof(GF_ELEMENT));
        memcpy(dec_ctx->message[pivot], message,  pktsize*sizeof(GF_ELEMENT));
        //printf("received-DoF %d new-DoF %d row_ops: %d\n", dec_ctx->DoF, pivot, rowop);
        dec_ctx->DoF += 1;
        if (!applying_precode)
            dofcount += 1;    // don't count dofs provided by the parity-check packets
    }
    return pivot;
}


// Apply the parity-check matrix to the decoding matrix
static int apply_parity_check_matrix(struct bats_decoder_ref *dec_ctx)
{
    static char fname[] = "apply_parity_check_matrix";
    int i, j, k;
    int num_of_new_DoF = 0;

    int pktsize = dec_ctx->param->pktsize;
    int numpp = dec_ctx->param->snum + dec_ctx->param->cnum;

    // 1, Copy parity-check vectors to the nonzero rows of the decoding matrix
    GF_ELEMENT *ces = malloc(numpp*sizeof(GF_ELEMENT));
    GF_ELEMENT *msg = malloc(pktsize*sizeof(GF_ELEMENT));
    int p = 0;          // index pointer to the parity-check vector that is to be copyed
    for (int p=0; p<dec_ctx->param->cnum; p++) {
        memset(ces, 0, numpp*sizeof(GF_ELEMENT));
        memset(msg, 0, pktsize*sizeof(GF_ELEMENT));
        /* Set the coding vector according to parity-check bits */
        NBR_node *varnode = dec_ctx->graph->l_nbrs_of_r[p]->first;
        while (varnode != NULL) {
            ces[varnode->data] = varnode->ce;
            varnode = varnode->next;
        }
        ces[dec_ctx->param->snum+p] = 1;
        int pivot = process_vector(dec_ctx, ces, msg);
    }
    free(ces);
    free(msg);

    /* Count available innovative rows */
    int missing_DoF = 0;
    for (i=0; i<numpp; i++) {
        if (dec_ctx->row[i] == NULL)
            missing_DoF += 1;
        else if (dec_ctx->row[i]->elem[0] ==0) {
            printf("%s: row[%d]->elem[0] is 0\n", fname, i);
        }
    }
    return missing_DoF;
}


static void back_substitution(struct bats_decoder_ref *dec_ctx)
{
    int pktsize = dec_ctx->param->pktsize;
    int numpp = dec_ctx->param->snum + dec_ctx->param->cnum;
    int i, j;
    int len;
    GF_ELEMENT quotient;
    for (i=numpp-1; i>=0; i--) {
        /* eliminate all nonzeros above diagonal elements from right to left*/
        for (j=0; j<i; j++) {
            len = dec_ctx->row[j]->len;
            if (j+len <= i || dec_ctx->row[j]->elem[i-j] == 0)
                continue;
            quotient = galois_divide(dec_ctx->row[j]->elem[i-j], dec_ctx->row[i]->elem[0]);
            galois_multiply_add_region(dec_ctx->message[j], dec_ctx->message[i], quotient, pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->row[j]->elem[i-j] = 0;
        }
        /* convert diagonal to 1*/
        if (dec_ctx->row[i]->elem[0] != 1) {
            galois_multiply_region(dec_ctx->message[i], galois_divide(1, dec_ctx->row[i]->elem[0]), pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->row[i]->elem[0] = 1;
        }
        /* save decoded packet */
        dec_ctx->pp[i] = calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dec_ctx->pp[i], dec_ctx->message[i], pktsize*sizeof(GF_ELEMENT));
    }
    dec_ctx->finished = 1;
}
