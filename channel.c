#include <stdio.h>
#include <stdlib.h>
#include "channel.h"


struct channel *create_channel(int delay, double pe)
{
    struct channel *ch = malloc(sizeof(struct channel));
    ch->delay = delay;
    ch->pe    = pe;
    ch->queue = calloc(delay+1, sizeof(void*));
    return ch;
}

void free_channel(struct channel *chnl)
{
    free(chnl->queue);
    free(chnl);
    chnl = NULL;
}

// Return whether the packet will be erased. If erased, it's the
// caller's responsibility to free the packet, as we have no idea
// about the packet's implementation (we always treat it as void)
int send_to_channel(struct channel *chnl, void *packet, int t)
{
    int queue_len = chnl->delay + 1;
    int lost = 0;
    int pos = (t+chnl->delay) % queue_len;      // where to put the packet in the queue 
    if (chnl->queue[pos] != NULL) {
        free(chnl->queue[pos]);                 // Warning: there might be memory leak here, but don't care it
    }
    if (rand() % 10000 <= chnl->pe * 10000) {
        lost = 1;
        chnl->queue[pos] = NULL;
        // printf("Element sent at time %d is lost\n", t);
    } else {
        chnl->queue[pos] = packet;
    }    
    return lost;
}

void *recv_from_channel(struct channel *chnl, int t)
{
    if (t < chnl->delay) {
        return NULL;
    }
    int queue_len = chnl->delay + 1;
    int pos = t % queue_len;      // where to pick the packet in the queue 
    void *out = chnl->queue[pos];
    chnl->queue[pos] = NULL;
    return out;
}

// Modify channel characteristics
// Important to realloc memory for queue, and arrange stored elements according to time
void modify_channel(struct channel *chnl, int t, int newdelay, double newpe)
{
    if (newdelay != chnl->delay) {
        void **newqueue = calloc(newdelay+1, sizeof(void *));
        int newqueue_len = newdelay + 1;
        int olddelay = chnl->delay;
        int oldqueue_len = chnl->delay + 1;
        if (newdelay < olddelay) {
            for (int i=0; i<=newdelay; i++) {
                int newpos = (t - i + newdelay) % newqueue_len;
                int oldpos = (t - i + olddelay) % oldqueue_len;
                newqueue[newpos] = chnl->queue[oldpos]; 
            }
        } else {
            for (int i=0; i<=olddelay; i++) {
                int newpos = (t - i + newdelay) % newqueue_len;
                int oldpos = (t - i + olddelay) % oldqueue_len;
                newqueue[newpos] = chnl->queue[oldpos]; 
            }

        }
        void **tmp = chnl->queue;
        chnl->queue = newqueue;
        chnl->delay = newdelay;
        free(tmp);
    }
    chnl->pe = newpe;
}