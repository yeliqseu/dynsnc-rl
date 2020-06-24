#include <stdlib.h>
// node of singly linked list
typedef struct node {
    int data;
    unsigned char ce;     // coefficient associated with node's data (for bipartite graph)
    struct node *next;
} NBR_node;

typedef struct node_list {
    struct node *first;
    struct node *last;
} NBR_nodes;

// Bipartitle graph
typedef struct bipartite_graph {
    int         nleft;
    int         nright;
    int         binaryce;       // Whether coefficients of edges are 1 or higher order
    NBR_nodes **l_nbrs_of_r;    // left side neighbours of right
    NBR_nodes **r_nbrs_of_l;    // right side neighbours of left
} BP_graph;

/* bipartite.c */
int create_bipartite_graph(BP_graph *graph, int nleft, int nright);
void free_bipartite_graph(BP_graph *graph);
