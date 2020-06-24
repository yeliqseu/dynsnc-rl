#ifndef BATS_H
#define BATS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bipartite.h"

#define ALIGN(a, b) ((a) % (b) == 0 ? (a)/(b) : (a)/(b) + 1)
typedef unsigned char GF_ELEMENT;

typedef struct bats_parameter {
    int         datasize;           // data size
    int         snum;               // number of source packets
    int         cnum;               // number of parity-check packets
    int         pktsize;            // packet content size
    int         seed;               // RNG seed
} BATSparam;

typedef struct bats_batch {
    int         batchid;            // batch id
    int         degree;             // number of packets in current batch (i.e., degree)
    int         bts;                // batch transmission size
    int         sent;               // number of packets sent from batch
    int         *pktid;             // packet id of the batch
} BATSbatch;

typedef struct bats_encoder_context {
    BATSparam   *param;
    BP_graph    *graph;
    int         batnum;             // number of batches generated so far
    BATSbatch   *currbat;           // current batch
    GF_ELEMENT  **pp;               // pointers to precoded source packets
} BATSencoder;

// Coded packet of bats code
typedef struct bats_packet {
    int         batchid;
    int         degree;             // number of packets encoded in this batch
    int         bts;                // batch transmission size
    int         *pktid;             // packet id of the packets
    GF_ELEMENT  *coes;              // coding coefficients
    GF_ELEMENT  *syms;              // coded packet content
} BATSpacket;

// Encoder
BATSencoder *bats_create_encoder(unsigned char *buf, BATSparam *param);
BATSbatch *bats_start_new_batch(BATSencoder *ctx, int batchid, int degree, int bts);
BATSpacket *bats_encode_packet(BATSencoder *ctx);
void bats_encode_packet_im(BATSencoder *ctx, BATSpacket *pkt);
BATSpacket *bats_duplicate_packet(BATSencoder *ctx, BATSpacket *pkt);
BATSpacket *bats_alloc_batch_packet(BATSencoder *ctx);
void bats_free_batch(BATSbatch *batch);
void bats_free_packet(BATSpacket *pkt);
void bats_free_encoder(BATSencoder *ctx);


// Recoding

// BATS buffer of fixed size
// Received packets from upstream are buffered in a first-in-first-out manner. 
// If the buffer is full, the oldest buffered packet would be discarded to store 
// a new received packet. We assume that the batch size is much larger than the 
// buffer size. Therefore, the buffered packets may belong to at most two different 
// batches, which happens when the buffer starts to receive the first several packets 
// belonging to a new batch. When there are two different batches in the buffer, 
// recoded packets are still generated from the older batch until all of its packets 
// are discarded to accommodate received packets of the new batch. The current batch 
// from which recoded packets are generated is referred to as the \textit{sending} batch 
// of the buffer, and the \textit{receiving} batch is the batch the latest received 
// packet belongs to. Clearly, the sending and receiving batches are the same if there 
// are only one batch in the buffer.
typedef struct bats_buffer {
    BATSparam   *param;             // pointer to the parameter of the BATS code 
    BATSpacket  **srbuf;
    int         bufsize;            // size of buffer
    int         sbatchid;           // current sending batch
    int         currbts;            // bts of the current sending batch
    int         s_first;            // start pos index of sending buffer
    int         r_last;              // end pos index of receiving buffer
                                    // if ((r_end+1) % bufsize == s_start), discard old pkt
} BATSbuffer;

BATSbuffer *bats_create_buffer(BATSparam *param, int bufsize);
void bats_buffer_packet(BATSbuffer *buf, BATSpacket *pkt);
BATSpacket *bats_recode_packet(BATSbuffer *buf);
void visualize_buffer(BATSbuffer *buf);
void bats_free_buffer(BATSbuffer *buf);

// Reference decoder

// Row vector of a matrix
struct row_vector
{
    int len;            // length of the row
    GF_ELEMENT *elem;   // elements of the row
} ROW;

// Reference decoder context
struct bats_decoder_ref {
    BATSparam           *param;             // BATS code parameter
    BP_graph            *graph;             // bipartite graph of precode

    int                 received;           // received coded packets
    int                 overhead;           // received coded packets
    int                 DoF;                // innovative coded packets
    int                 *seen;              // 0/1 array to indicate whether a packet has been included in at least one received batch
    int                 covered;            // covered packets in the received batches
    int                 de_precode;         // whether applied parity-check vectors?
    int                 finished;           // whether finished decoding
    long long           operations;         // finite field operations 
    struct row_vector   **row;              // rows of decoding matrix
    GF_ELEMENT          **message;          // rows of message symbols
    GF_ELEMENT          **pp;               // recovered packets

    // A small matrix storing vectors of the receiving batch. Received vectors are processed
    // against the previously vectors of the same batch, which renders the vector sparser.
    int                 currbid;
    int                 currpnum;
    struct row_vector **batch_row;
    GF_ELEMENT        **batch_msg;
};

struct bats_decoder_ref *bats_create_decoder_ref(BATSparam *param);
void bats_free_decoder_ref(struct bats_decoder_ref *decoder);
int bats_process_packet_ref(struct bats_decoder_ref *dec_ctx, BATSpacket *pkt);

// Reinforcement learning functions
int derive_e_greedy_action_SGD(double r_ratio, int isgreedy);
double calculate_action_value_q_estimate(double r_ratio, int act_id, int tiles_array[]);
int simulate_with_current_weights(int episode, int nsim, int snum, int cnum, int nhop, int bufsize);
#endif
