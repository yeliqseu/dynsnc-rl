#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "galois.h"
#include "bats.h"

extern void init_genrand(unsigned long s);
extern unsigned long genrand_int32(void);
static void bats_precoding(BATSencoder *ctx);
static int compare_int(const void *elem1, const void *elem2);
static void get_random_unique_numbers(int ids[], int n, int ub);

BATSencoder *bats_create_encoder(unsigned char *buf, BATSparam *param)
{
    static char fname[] = "bats_create_encoder";
    // Allocate encoder context
    BATSencoder *ctx;
    if ( (ctx = calloc(1, sizeof(BATSencoder))) == NULL ) {
        fprintf(stderr, "%s: calloc encctx\n", fname);
        return NULL;
    }
    ctx->param = param;
    // init mt19937 RNG
    init_genrand(param->seed);
    // calculate number of source packets after padding 0 (in case)
    param->snum = ALIGN(param->datasize, param->pktsize);
    // create bipartite graph
    if (param->cnum != 0) {
        if ( (ctx->graph = malloc(sizeof(BP_graph))) == NULL ) {
            fprintf(stderr, "%s: malloc BP_graph\n", fname);
            return NULL;
        }
        if (create_bipartite_graph(ctx->graph, param->snum, param->cnum) < 0)
            return NULL;
    }
    ctx->batnum = 0;
    ctx->currbat = NULL;

    constructField();   // Construct Galois Field

    // load source packets and perform precoding
    if ((ctx->pp = calloc(param->snum+param->cnum, sizeof(GF_ELEMENT*))) == NULL) {
        fprintf(stderr, "%s: calloc ctx->pp\n", fname);
        //bats_free_encoder(ctx);
        return NULL;
    }

    if (buf != NULL) {
        int alread = 0;
        int i;
        // Load source packets
        for (i=0; i<param->snum; i++) {
            ctx->pp[i] = calloc(param->pktsize, sizeof(GF_ELEMENT));
            int toread = (alread+param->pktsize) <= param->datasize ? param->pktsize : param->datasize-alread;
            memcpy(ctx->pp[i], buf+alread, toread*sizeof(GF_ELEMENT));
            alread += toread;
        }
        // Allocate parity-check packet space
        for (i=0; i<param->cnum; i++)
            ctx->pp[param->snum+i] = calloc(param->pktsize, sizeof(GF_ELEMENT));
        bats_precoding(ctx);
    }


    return ctx;
}


static void bats_precoding(BATSencoder *ctx)
{
    static char fname[] = "perform_precoding";

    int i, j;
    for (i=0; i<ctx->param->cnum; i++) {
        // Encoding check packet according to the LDPC graph
        NBR_node *nb = ctx->graph->l_nbrs_of_r[i]->first;
        while(nb != NULL) {
            int sid = nb->data;  // index of source packet
            // XOR information content
            galois_multiply_add_region(ctx->pp[i+ctx->param->snum], ctx->pp[sid], nb->ce, ctx->param->pktsize);
            // move to next possible neighbour node of current check
            nb = nb->next;
        }
    }
}

// Generate a new batch from soruce/intermediate packets
BATSbatch *bats_start_new_batch(BATSencoder *ctx, int batchid, int degree, int bts)
{
    static char fname[] = "bats_new_batch";
    if (ctx->currbat != NULL) {
        //printf("Encoder has a current batch, let me free it first.\n");
        bats_free_batch(ctx->currbat);
    }

    BATSbatch *batch = malloc(sizeof(BATSbatch));
    batch->batchid = batchid;
    batch->degree = degree;
    batch->bts    = bts;
    batch->sent = 0;
    batch->pktid = malloc(sizeof(int)*batch->degree);
    // uniformatly randomly draw packets from source/intermediate packets using Fisher-Yates algo.
    get_random_unique_numbers(batch->pktid, batch->degree, ctx->param->snum+ctx->param->cnum);
    ctx->currbat = batch;
    ctx->batnum += 1;
    return batch;
}

// generate a number of n<ub unique random numbers within the range of [0, ub-1] 
// using Fisher-Yates shuffle method
static void get_random_unique_numbers(int ids[], int n, int ub)
{
	int init_array[ub];
	int i, j;
	for (i=0; i<ub; i++)
		init_array[i] = i;

	// randomly shuffle the init_array
	for (i=ub-1; i>=1; i--) {
		int rd = genrand_int32() % (i+1);
		//int rd = gsl_rng_uniform_int(r, i+1);
		int tmp = init_array[rd];
		init_array[rd] = init_array[i];
		init_array[i] = tmp;
	}

    // sort the obtained unique random numbers so that coding coefficients corresponding
    // to packets are stored in the ascending order (to simplify decoder implementation)
    qsort(init_array, n, sizeof(int), compare_int);
    memcpy(ids, init_array, n*sizeof(int));
	//for (j=0; j<n; j++)
	//	ids[j] = init_array[j];
}

static int compare_int(const void *elem1, const void *elem2)
{
    int a = * ((int *) elem1);
    int b = * ((int *) elem2);
    if (a < b)
        return -1;
    if (a > b)
        return 1;
    return 0;
}

// Encode a packet
BATSpacket *bats_encode_packet(BATSencoder *ctx)
{
    // allocate memory
    BATSpacket *newpkt = bats_alloc_batch_packet(ctx);
    bats_encode_packet_im(ctx, newpkt);
    return newpkt;
}

// Encode a packet (already in memory) from current batch
void bats_encode_packet_im(BATSencoder *ctx, BATSpacket *pkt)
{
    pkt->batchid = ctx->currbat->batchid;
    int n = ctx->currbat->degree;
    pkt->degree  = n;
    memcpy(pkt->pktid, ctx->currbat->pktid, sizeof(int)*n);
    pkt->bts = ctx->currbat->bts;
    // start encoding
    GF_ELEMENT co;
    memset(pkt->syms, 0, sizeof(GF_ELEMENT)*ctx->param->pktsize);
    for (int i=0; i<n; i++) {
        int id = ctx->currbat->pktid[i];
        co = (GF_ELEMENT) genrand_int32() % (1 << 8);     // Randomly generated coding coefficient
        pkt->coes[i] = co;
        galois_multiply_add_region(pkt->syms, ctx->pp[id], co, ctx->param->pktsize);
    }
    ctx->currbat->sent += 1;
}

// Allocate empty BATSpacket
BATSpacket *bats_alloc_batch_packet(BATSencoder *ctx)
{
    BATSpacket *pkt = calloc(1, sizeof(BATSpacket));
    if (pkt == NULL)
        return NULL;
    pkt->batchid = ctx->currbat->batchid;

    int degree = ctx->currbat->degree;
    pkt->degree = degree;
    pkt->pktid = calloc(degree, sizeof(int));
    if (pkt->pktid == NULL)
        goto AllocErr;
    memcpy(pkt->pktid, ctx->currbat->pktid, sizeof(int)*degree);
    pkt->coes = calloc(degree, sizeof(GF_ELEMENT));
    if (pkt->coes == NULL)
        goto AllocErr;
    pkt->syms = calloc(ctx->param->pktsize, sizeof(GF_ELEMENT));
    if (pkt->syms == NULL)
        goto AllocErr;

    return pkt;

AllocErr:
    bats_free_packet(pkt);
    return NULL;
}

BATSpacket *bats_duplicate_packet(BATSencoder *ctx, BATSpacket *pkt)
{
    BATSpacket *dup_pkt = bats_alloc_batch_packet(ctx);
    dup_pkt->batchid = pkt->batchid;
    dup_pkt->degree = pkt->degree;
    memcpy(dup_pkt->pktid, pkt->pktid, sizeof(int)*pkt->degree);
    memcpy(dup_pkt->coes, pkt->coes, sizeof(GF_ELEMENT)*pkt->degree);
    memcpy(dup_pkt->syms, pkt->syms, sizeof(GF_ELEMENT)*ctx->param->pktsize);
    return dup_pkt;
}


void bats_free_batch(BATSbatch *batch)
{
    if (batch->pktid != NULL)
        free(batch->pktid);
    free(batch);
    batch = NULL;
    return;
}

void bats_free_encoder(BATSencoder *ctx)
{
    free_bipartite_graph(ctx->graph);
    ctx->graph = NULL;
    if (ctx->currbat != NULL) {
        bats_free_batch(ctx->currbat);
        ctx->currbat = NULL;
    }
    for (int i=0; i<ctx->param->snum+ctx->param->cnum; i++) {
        free(ctx->pp[i]);
        ctx->pp[i] = NULL;
    }
    free(ctx->pp);
    // free(ctx->param);
}

void bats_free_packet(BATSpacket *pkt)
{
    if (pkt == NULL)
        return;
    if (pkt->pktid != NULL)
        free(pkt->pktid);
    if (pkt->coes != NULL)
        free(pkt->coes);
    if (pkt->syms != NULL)
        free(pkt->syms);
    free(pkt);
    return;
}
