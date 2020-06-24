#include "bats.h"
#include "galois.h"


extern unsigned long genrand_int32(void);
static int compare_int(const void *elem1, const void *elem2);
static void get_random_unique_numbers(int ids[], int n, int ub);
// Synchronized recoding: 
// 1st and 2nd hop always send packets of the same batch.

static int s_neq_r = 0;  // indicate whether the sending and receiving batch do not match.
                     // i.e., whether there are more than 1 batches buffered


static int systematic = 0;       // in the systematic phase
static int newsys = -1;          // an uncoded packet is received, can be forwarded (the value is the storage index in the buffer)
static int recoded_sys = 9999;   // the "batchid" for the coded packets generated from the buffered uncoded packets

extern int s_count;  // count num of pkts sent from the current sending batch
// extern int currbatch;

void visualize_buffer(BATSbuffer *buf);

BATSbuffer *bats_create_buffer(BATSparam *param, int bufsize)
{
    static char fname[] = "bats_create_buffer";

    char *syst = getenv("BATS_SYSTEMATIC");
    if (syst != NULL)
        systematic = atoi(syst);

    BATSbuffer *buf = malloc(sizeof(BATSbuffer));

    buf->param = param;
    buf->srbuf = calloc(bufsize, sizeof(BATSpacket *));  // allocate buffer packet pointers
    buf->bufsize = bufsize;
    buf->sbatchid = -1;  // empty buffer
    buf->currbts = -1;
    buf->s_first = -1;
    buf->r_last = -1;
    return buf;
}

void bats_buffer_packet(BATSbuffer *buf, BATSpacket *pkt)
{
    int pos = -1;   // Pos index where new packet is stored

    // a new batch has started, flush the buffer
    if (pkt->batchid != buf->sbatchid) {
        printf("new batch %d is seen while the current recoding batch is %d\n", pkt->batchid, buf->sbatchid);
        buf->sbatchid = pkt->batchid;
        buf->currbts = pkt->bts;
        // overwrite previous bufferred packets
        pos = 0;
        buf->r_last  = pos;
        buf->s_first = pos;
        if (buf->srbuf[pos] != NULL)
            bats_free_packet(buf->srbuf[pos]);
        buf->srbuf[pos] = pkt;
        return;
    }

    if (((buf->r_last+1) % buf->bufsize) == buf->s_first) {
        // buffer is full, replace the oldest packet with the new arrived one
        int new_sfirst = (buf->s_first+1) % buf->bufsize;

        //visualize_buffer(buf);
        
        pos = buf->s_first;
        bats_free_packet(buf->srbuf[pos]);
        buf->srbuf[pos] = NULL;
        buf->srbuf[pos] = pkt;
        buf->r_last = pos;
        buf->s_first = new_sfirst;   // update the head coded packet of the batch
        return;
    } else {
        pos = (buf->r_last + 1) % buf->bufsize;
        // just double check
        if (buf->srbuf[pos] != NULL) {
            bats_free_packet(buf->srbuf[pos]);
            buf->srbuf[pos] = NULL;
        }
        buf->srbuf[pos] = pkt;
        buf->r_last = pos;
        return;
    }

    return ;
}

BATSpacket *bats_recode_packet(BATSbuffer *buf)
{
    int i;

    if (buf->sbatchid < 0) {
        //printf("Buffer has no batch buffered yet\n");
        return NULL;
    }

    int s_pos = buf->s_first; // first packet of the sending batch
    int pos;
    

    BATSpacket *pkt = malloc(sizeof(BATSpacket));
    pkt->batchid = buf->sbatchid;
    //printf("recoded batch: %d\n", pkt->batchid);
    pkt->degree = buf->srbuf[s_pos]->degree;
    pkt->bts = buf->srbuf[s_pos]->bts;
    pkt->pktid = malloc(sizeof(int)*pkt->degree);
    memcpy(pkt->pktid, buf->srbuf[s_pos]->pktid, sizeof(int)*pkt->degree);

    pkt->coes = calloc(pkt->degree, sizeof(GF_ELEMENT));
    pkt->syms = calloc(buf->param->pktsize, sizeof(GF_ELEMENT));
    GF_ELEMENT co =0;    
    for (i=0; i<buf->bufsize; i++) {
        pos = (s_pos + i) % buf->bufsize;
        if (buf->srbuf[pos] == NULL || buf->srbuf[pos]->batchid != buf->sbatchid)
            break;      // packets belonging to the same batch must be stored adjacently.
        co = genrand_int32() % (1<<8);
        galois_multiply_add_region(pkt->coes, buf->srbuf[pos]->coes, co, pkt->degree);
        galois_multiply_add_region(pkt->syms, buf->srbuf[pos]->syms, co, buf->param->pktsize);
    }
    // s_count += 1;
    return pkt;
}

void visualize_buffer(BATSbuffer *buf)
{
    printf("buffer size: %d sbatchid: %d s_first: %d r_last: %d s_count: %d s_neq_r: %d\n", buf->bufsize, buf->sbatchid, buf->s_first, buf->r_last, s_count, s_neq_r);
    for (int i=0; i<buf->bufsize; i++) {
        if (buf->srbuf[i] != NULL) {
            printf("%d\t", buf->srbuf[i]->batchid);
        } else {
            printf("-1\t");
        }
    }
    printf("\n");
}

void bats_free_buffer(BATSbuffer *buf)
{
    for (int i=0; i<buf->bufsize; i++) {
        if (buf->srbuf[i] != NULL) {
            bats_free_packet(buf->srbuf[i]);
            buf->srbuf[i] = NULL;
        }
    }
    free(buf->srbuf);
    buf->srbuf = NULL;
    free(buf);
    buf = NULL;
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
