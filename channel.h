#include <stdio.h>

/*
 * Channel with delay and erasure
 * Implemented using an array to simulate a delayed queue
 */

typedef struct channel {
    void    **queue;
    int     delay;
    double  pe;
} Channel;

struct channel *create_channel(int delay, double pe);
void free_channel(struct channel *chnl);
int send_to_channel(struct channel *chnl, void *packet, int t);
void *recv_from_channel(struct channel *chnl, int t);

// The inflight packets of the channel is implemented as an array of pointers (queue) of length equal
// to delay+1. If the new delay is smaller than the old, the queue would be shrinked. The caller of the
// channel has to ensure that recv_from_channle() is called for (olddelay-newdelay) times to receive the
// packets in the shrinked part of the old queue.
void modify_channel(struct channel *chnl, int t, int delay, double pe);

