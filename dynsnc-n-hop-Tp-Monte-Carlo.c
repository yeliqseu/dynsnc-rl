// Impliment on-policy Monte Carlo control for dynamic SNC
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

double alpha;
double beta;
double gamma;
double epsilon;
double **Qtable;
double **Rtable;
double **piPolicy;
double **visits;
int    nstate;
int    naction = 3;
int    action[3][2] = {
                        { 16, 4},
                        { 16, 8},
                        { 16, 16},
                    };

double pe;

static int derive_action_id_from_policy(int r);

char usage[] = "Simulate n-hop lossy line networks\n\
                \n\
                S ---(1-pe_1)---> V1 ---(1-pe_2)---> ... ---(1-pe_n)---> D\n\
                \n\
                usage: ./programName nslots alpha beta gamma epsilon snum cnum pe nhop \n\
                nslots   - number of simulation time slots\n\
                alpha    - learning step size\n\
                beta     - reward parameter\n\
                gamma    - discounting rate\n\
                epsilon  - epsilon-greedy parameter\n\
                snum     - number of source packets\n\
                cnum     - number of parity-check packets of precode\n\
                pe       - initial packet loss probability on each hop (equal)\n\
                nhop     - number of hops (integer)\n\
                Tp       - propagation delay on each hop (equal)\n\
                lambda   - lambda of Q(lambda)\n";
int main(int argc, char *argv[])
{
    if (argc != 12) {
        printf("%s\n", usage);
        exit(1);
    }

    int i, j;
    // Q-learning parameters
    int nslots = atoi(argv[1]);
    alpha = atof(argv[2]);
    beta = atof(argv[3]);
    gamma = atof(argv[4]);
    epsilon = atof(argv[5]);
    // Network and coding parameters
    int snum = atoi(argv[6]);
    int cnum = atoi(argv[7]);
    pe = atof(argv[8]);
    int nhop = atoi(argv[9]);
    int Tp   = atoi(argv[10]);
    double lambda = atof(argv[11]);
    int bufsize = snum+cnum;
    int pktsize = 256;
    int datasize = snum * pktsize;
    
    // Initialize Q table
    // Allocate memory, and then try to load existing Q-table and learn based on it
    int K = snum + cnum;
    //nstate = (K + 1) * (K + 2) / 2;
    nstate = K + 1;
    Qtable = calloc(nstate, sizeof(double*));
    Rtable = calloc(nstate, sizeof(double*));
    piPolicy = calloc(nstate, sizeof(double*));
    for (i=0; i<nstate; i++) {
        Qtable[i] = calloc(naction, sizeof(double));
        Rtable[i] = calloc(naction, sizeof(double));
        piPolicy[i] = calloc(naction, sizeof(double)); 
    }
    // Initialize Q-table arbitrarily
    // And initialize an arbitrary e-soft policy (epsilon for others, and (1 - |A|*epsilon) for one arbitrary state
    for (i=0; i<nstate; i++) {
        for (int j=0; j<naction; j++) {
            Qtable[i][j] =  - (double) rand() / RAND_MAX;
            piPolicy[i][j] = epsilon;
        }
        piPolicy[i][rand() % naction] = 1 - (naction - 1) * epsilon;
    }
    // Allocate/load table for counting the number of visits to each entry of the Q-table
    visits = calloc(nstate, sizeof(double*));
    for (i=0; i<nstate; i++) {
        visits[i] = calloc(naction, sizeof(double));
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec

    // generate random source data
    unsigned char *databuf = malloc(datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, databuf, datasize);
    close(rnd);
    int t = 0;
    int episode = 0;
    int newgen = 1;         // indicate a generation of transmissions
    while (t < nslots) {
        if (newgen) {
            episode += 1;
            printf("Learning from episode %d...\n", episode);
            newgen = 0;
            //epsilon = epsilon * 0.99;
        }

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


        // create channels
        struct channel **chnl = calloc(nhop, sizeof(struct chnl*));
        double fb_succ = 1.0;
        for (i=0; i<nhop; i++) {
            chnl[i] = create_channel(Tp, pe);
            fb_succ *= (1-pe);
        }
        // create feedback channel (note multihop relay) 
        struct channel *fdbk_chnl;
        char *perfect_fb = getenv("PERFECT_FB");
        if (perfect_fb != NULL && strcmp(perfect_fb, "TRUE") == 0) {
            fdbk_chnl = create_channel(0, 0);
        } else {
            fdbk_chnl = create_channel(Tp*nhop, 1-fb_succ);
        }
        
        // arrays to store information for processing after each episode
        int *state_seq = calloc(10000, sizeof(int));
        int *action_seq = calloc(10000, sizeof(int));
        int *reward_seq = calloc(10000, sizeof(int));
        int act_seq = 0;

        BATSpacket *pkt;
        int nuse = 0;                   // number of network uses
        int r_curr = 0;
        int r_next = 0;
        //int act_id = derive_e_greedy_action(r_curr, Qtable, naction);
        int act_id = derive_action_id_from_policy(r_curr);  // derive from the e-soft policy
        int newDeg = action[act_id][0];
        int newBts = action[act_id][1];
        // create new batch
        printf("Current (r, c) = ( %d %d ). Action: deg= %d bts= %d\n", 0, 0, newDeg, newBts);
        bats_start_new_batch(encoder, currbatch, newDeg, newBts);

        while (!decoder->finished) {
            // change channel paramter if an environment argument is set
            if (t== nslots/3 || t == 2 * nslots/3) {
                char *changing = getenv("NONSTATIONARY");
                if (changing != NULL && strcmp(changing, "TRUE") == 0) {
                    if (t == nslots/3) {
                            pe = pe * 2;
                            Tp = Tp * 2;
                    }
                    if (t == 2 * nslots /3) {
                            pe = pe / 4;
                            Tp = Tp / 4;
                    }
                    // modify channel parameters
                    printf("At time %d, channel paramters are changed\n", t);
                    fb_succ = 1.0;
                    for (i=0; i<nhop; i++) {
                        modify_channel(chnl[i], t, Tp, pe);
                        fb_succ *= (1-pe);
                    }
                    char *perfect_fb = getenv("PERFECT_FB");
                    if (perfect_fb != NULL && strcmp(perfect_fb, "TRUE") == 0) {
                        modify_channel(fdbk_chnl, t, 0, 0);
                    } else {
                        modify_channel(fdbk_chnl, t, Tp *2, 1-fb_succ);
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
            // decoder feedback to sender its currently received DoF
            int *r_feedback = calloc(1, sizeof(int));
            r_feedback[0] = decoder->DoF;
            printf("Decoder send feedback r = %d at time %d\n", decoder->DoF, nuse);
            int fblost = send_to_channel(fdbk_chnl, r_feedback, nuse);
            if (fblost) {
                free(r_feedback);
                r_feedback= NULL;
            }
            // sender receive r feedback
            int *recv_r = (int *) recv_from_channel(fdbk_chnl, nuse);
            if (recv_r != NULL) {
                r_next = recv_r[0];         // update r knowledge at the sender
                printf("Received feedback r = %d from decoder at time %d\n", r_next, nuse);
                free(recv_r);
                recv_r = NULL;
            }

            // check whether it's time to change a batch
            if (batchsent >= encoder->currbat->bts || decoder->finished) {
                // will start a new batch, summarize last batch statistics
                // printf("Received %d packets from batch %d (BTS: %d) and %d are innovative [sent from buffer: %d]\n", batchcount, currbatch, encoder->currbat->bts, dofcount, s_count);
                // printf("Batch %d is done, start next batch...\n", currbatch);
                // get (r, c) feedback from decoder, i.e., the observed new state
                int oldDeg = encoder->currbat->degree;
                int oldBts = encoder->currbat->bts;
                double reward = 0.0;
                if (decoder->finished) {
                    reward = - batchsent;
                } else {
                    reward = - oldBts;
                }
                // Record state, action, and reward sequence
                state_seq[act_seq] = r_curr;
                action_seq[act_seq] = act_id;
                reward_seq[act_seq] = reward;
                
                printf(" r_prev= %d , batch %d action: deg= %d bts= %d r_new= %d , reward: %.4f next action: deg= %d bts=%d\n", 
                         r_curr, currbatch, oldDeg, oldBts, r_next, reward, newDeg, newBts);
                // choose NEW action following policy pi
                r_curr = r_next;
                act_id = derive_action_id_from_policy(r_curr);
                newDeg = action[act_id][0];
                newBts = action[act_id][1];
                act_seq += 1;

                // Start next batch if the episode is not yet complete
                if (!decoder->finished) {
                    currbatch++;
                    bats_start_new_batch(encoder, currbatch, newDeg, newBts);
                }
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

        printf("time: %d snum: %d cnum: %d pktsize: %d nbatch: %d overhead: %.4f ops: %.6f network-uses: %d ( episode: %d )\n", 
               t, param.snum, param.cnum, param.pktsize, encoder->batnum, 
               (double) decoder->overhead/decoder->param->snum, (double) decoder->operations/decoder->param->snum/decoder->param->pktsize, nuse, episode);

        // process the previous episode, backward loop
        double G = 0.0;
        int s_prev = -1;
        int a_prev = -1;
        for (int i=act_seq-1; i--; i>=0) {
            G = gamma * G + reward_seq[i];
            int s = state_seq[i];
            int a = action_seq[i];
            if (s == s_prev && a == a_prev) {
                continue;
                // Note: since state of our problem is non-decreasing, we only need to compare with the previous
                // state and action to ensure first-visit in the current episode.
            }
            Rtable[s][a] += G;
            visits[s][a] += 1;
            Qtable[s][a] = Rtable[s][a] / visits[s][a];
            int opt_actid = derive_optimal_action(s, Qtable, naction);
            for (int j=0; j<naction; j++) {
                if (j == opt_actid) {
                    piPolicy[s][j] = 1 - epsilon + epsilon / naction;
                } else {
                    piPolicy[s][j] = epsilon / naction;
                }
            }
            s_prev = s;
            a_prev = a;
        }
        free(state_seq);
        free(reward_seq);
        free(action_seq);
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
        free_channel(fdbk_chnl);
        fdbk_chnl = NULL;
        // encoder
        batchsent = 0;              // number of packets sent from the current batch
        // recoder
        s_count = 0;    // number of recoded packets sent from the current batc
        // decoder
        batchcount = 0;   // number of received packets of the current receiving batch
        dofcount = 0;     // number of innovative packets contributed by the current receiving batch
        newgen = 1;     // to start a new generation (i.e., episode)
    }
    //save_table(Qfname, Qtable, nstate, naction);
    //save_table(vQfname, visits, nstate, naction);
    free(databuf);
    return 0;
}

static int derive_action_id_from_policy(int r)
{
    int x = rand() % 100000;
    int sum = 0;
    int state_id = r;
    int act_id = 0;
    for (int i=0; i<naction; i++) {
        sum += piPolicy[state_id][i] * 100000;
        if (x < sum) {
            act_id = i;
            break;
        }
    }
    return act_id;
}