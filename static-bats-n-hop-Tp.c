// Each hop is of a dynamically changing erasure probability
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>

#include "bats.h"
#include "channel.h"

int currbatch = 0;
// encoder
int batchsent = 0;              // number of packets sent from the current batch
// recoder
int s_count = 0;    // number of recoded packets sent from the current batc
// decoder
int batchcount = 0;   // number of received packets of the current receiving batch
int dofcount = 0;     // number of innovative packets contributed by the current receiving batch


char usage[] = "Simulate n-hop lossy line networks\n\
                \n\
                S ---(1-pe_1)---> V1 ---(1-pe_2)---> ... ---(1-pe_n)---> D\n\
                \n\
                usage: ./programName nslots snum cnum deg bts pe nhop Tp\n\
                nslots   - number of time slots for simulation\n\
                snum     - number of source packets\n\
                cnum     - number of parity-check packets of precode\n\
                deg      - degree of batch\n\
                bts      - bts of batch\n\
                pe       - initial packet loss probability on each hop (equal)\n\
                nhop     - number of hops (integer)\n\
                Tp       - propagation delay on each hop (equal)\n";
int main(int argc, char *argv[])
{
    if (argc != 9) {
        printf("%s\n", usage);
        exit(1);
    }

    int i;
    // Network and coding parameters
    int nslots = atoi(argv[1]);
    int snum   = atoi(argv[2]);
    int cnum   = atoi(argv[3]);
    int deg    = atoi(argv[4]);
    int bts    = atoi(argv[5]);
    double pe  = atof(argv[6]);
    int nhop   = atoi(argv[7]);
    int Tp     = atoi(argv[8]);
    int bufsize = snum+cnum;
    int pktsize = 256;
    int datasize = snum * pktsize;
    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec

    // generate random source data
    unsigned char *databuf = malloc(datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, databuf, datasize);
    close(rnd);
    int t = 0;
    while (t < nslots) {
        BATSparam param = { datasize,
                            0,              // no need to specify number of source packets
                            cnum,           // parity-check packets
                            pktsize,
                            0,              // seed for RNG
                        };
        // create encoder at source node
        BATSencoder *encoder = bats_create_encoder(databuf, &param);
        // create recoders at intermediate nodes
        BATSbuffer **buf = calloc(nhop-1, sizeof(BATSbuffer*));
        for (i=0; i<nhop-1; i++)
            buf[i] = bats_create_buffer(&param, bufsize);
        // create decoder at destination node
        struct bats_decoder_ref *decoder = bats_create_decoder_ref(&param);

        // create forward channels
        struct channel **chnl = calloc(nhop, sizeof(struct chnl*));
        for (i=0; i<nhop; i++) {
            chnl[i] = create_channel(Tp, pe);
        }

        BATSpacket *pkt;
        int nuse = 0;                   // number of network uses
        // create new batch
        bats_start_new_batch(encoder, currbatch, deg, bts);

        while (!decoder->finished) {
            // change channel paramter if an environment argument is set
            if (t== nslots/3 || t == 2 * nslots/3) {
                char *changing = getenv("NONSTATIONARY");
                if (changing != NULL && strcmp(changing, "TRUE") == 0) {
                    if (t == nslots/3) {
                            pe = pe * 2;
                            Tp = Tp * 2;
                            // modify channel parameters
                            printf("At time %d, channel paramters are changed\n", t);
                            for (i=0; i<nhop; i++) {
                                modify_channel(chnl[i], t, Tp, pe);
                            }
                    }
                    if (t == 2 * nslots /3) {
                            pe = pe / 4;
                            Tp = Tp / 4;
                            // modify channel parameters
                            printf("At time %d, channel paramters are changed\n", t);
                            for (i=0; i<nhop; i++) {
                                modify_channel(chnl[i], t, Tp, pe);
                            }
                    }
                }
            }
            // use each forward hop once
            for (int i=0; i<nhop; i++) {
                // first hop
                if (i == 0) {
                    pkt = bats_encode_packet(encoder);
                } else {
                    pkt = bats_recode_packet(buf[i-1]);
                    if (pkt == NULL) {
                        continue;
                    }
                }
                // send to channel of hop i
                int lost = send_to_channel(chnl[i], pkt, nuse);
                if (lost) {
                    bats_free_packet(pkt);
                    pkt = NULL;
                    printf("Packet sent on hop %d at time %d is lost\n", i, nuse);
                }
                // receive from channel of hop i
                BATSpacket *rpkt = (BATSpacket*) recv_from_channel(chnl[i], nuse);
                if (rpkt != NULL) {
                    if (i < nhop-1) {
                        bats_buffer_packet(buf[i], rpkt);        // next node is an intermediate node
                    } else {
                        bats_process_packet_ref(decoder, rpkt);  // next node is decoder
                    }
                } else {
                    printf("Receive NULL from channel of hop %d at time %d\n", i, nuse);
                }
            }

            // check whether it's time to change a batch
            if (batchsent >= encoder->currbat->bts) {
                currbatch++;
                bats_start_new_batch(encoder, currbatch, deg, bts);
                // reset batchsent
                batchsent = 0;
                batchcount = 0;
                dofcount = 0;
                s_count = 0;
            }
            t++;
            nuse++;
            batchsent++;
        }
        for (int i=0; i<param.snum; i++) {
            if (memcmp(encoder->pp[i], decoder->pp[i], param.pktsize) != 0) {
                fprintf(stderr, "recovered is NOT identical to original.\n");
                break;
            }
        }

        printf("bufsize: %d numhop: %d ", bufsize, nhop);
        printf("\n");

        printf("time: %d snum: %d cnum: %d pktsize: %d degree: %d bts: %d nbatch: %d overhead: %.4f ops: %.6f network-uses: %d \n", 
                t, param.snum, param.cnum, param.pktsize, deg, bts, encoder->batnum, 
                (double) decoder->overhead/decoder->param->snum, (double) decoder->operations/decoder->param->snum/decoder->param->pktsize, nuse);
        // free memory allocation
        // free encoder
        bats_free_encoder(encoder);
        encoder = NULL;
        // free recoder buffers
        for (i=0; i<nhop-1;i++) {
            if (buf[i] != NULL) {
                bats_free_buffer(buf[i]);
                buf[i] = NULL;
            }
        }
        free(buf);
        buf = NULL;
        // free decoder
        bats_free_decoder_ref(decoder);
        decoder = NULL;
        currbatch = 0;
        // free channel
        for (i=0; i<nhop; i++) {
            free_channel(chnl[i]);
            chnl[i] = NULL;
        }
        // encoder
        batchsent = 0;              // number of packets sent from the current batch
        // recoder
        s_count = 0;    // number of recoded packets sent from the current batc
        // decoder
        batchcount = 0;   // number of received packets of the current receiving batch
        dofcount = 0;     // number of innovative packets contributed by the current receiving batch
    }
    free(databuf);
    return 0;
}
