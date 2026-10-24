#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#define DIST_INF (LLONG_MAX/4)

typedef struct road{
    int t;
    int w;
}road;

typedef struct Graph{
    int n;
    int m;
    int *offset;
    road *edges;
}Graph;

typedef struct HeapItem{
    int vertex;
    long long key;
}HeapItem;

typedef struct BinaryHeap{
    HeapItem *items;
    int *position;
    int size;
    int capacity;
}BinaryHeap;

/* Array-backed 3-ary min-heap; position[v] stores v's array index. */
typedef struct TernaryHeap{
    HeapItem *items;
    int *position;
    int size;
    int capacity;
}TernaryHeap;

typedef struct FibNode{
    int vertex;
    long long key;
    int degree;
    int mark;
    int in_heap;
    struct FibNode *parent;
    struct FibNode *child;
    struct FibNode *left;
    struct FibNode *right;
}FibNode;

typedef struct FibHeap{
    FibNode *nodes;
    FibNode **root_buffer;
    FibNode **degree_table;
    int degree_capacity;
    int n;
    int size;
    FibNode *min;
}FibHeap;

typedef struct BinomialNode{
    int vertex;
    long long key;
    int degree;
    struct BinomialNode *parent;
    struct BinomialNode *child;
    struct BinomialNode *sibling;
}BinomialNode;

typedef struct BinomialHeap{
    BinomialNode *nodes;
    BinomialNode **node_of;
    BinomialNode *roots;
    int n;
    int size;
    int node_count;
}BinomialHeap;

typedef struct DijkstraWorkspace{
    long long *dist;
    unsigned char *state;
    int *touched;
    int touched_count;
    int n;
}DijkstraWorkspace;

typedef struct SearchStats{
    long long scanned_arcs;
    long long settled_vertices;
    long long pq_inserts;
    long long pq_extracts;
    long long pq_decrease_keys;
}SearchStats;

/* Internal heap counters. Comparisons count key/item comparator calls;
   binomial degree/list decisions are tracked separately. */
typedef struct HeapInternalStats{
    long long comparisons;
    long long upward_steps;
    long long downward_iterations;
    long long array_swaps;
    long long key_swaps;
    long long tree_links;
    long long cuts;
    long long cascade_steps;
    long long root_list_steps;
    long long children_promoted;
    long long max_queue_size;
    long long max_root_count;
    long long max_degree;
}HeapInternalStats;

/* Internal heap counters are compiled only into the counting build
   (-DCOUNT_HEAP_OPS). The timing build expands both macros to nothing,
   so no counter check or update sits in the measured path. */
#ifdef COUNT_HEAP_OPS
static HeapInternalStats *active_heap_metrics=NULL;

#define METRIC_ADD(field,value) do{ \
    if(active_heap_metrics!=NULL){ \
        active_heap_metrics->field+=(value); \
    } \
}while(0)

#define METRIC_MAX(field,value) do{ \
    if(active_heap_metrics!=NULL&& \
       active_heap_metrics->field<(long long)(value)){ \
        active_heap_metrics->field=(long long)(value); \
    } \
}while(0)
#else
#define METRIC_ADD(field,value) do{}while(0)
#define METRIC_MAX(field,value) do{}while(0)
#endif

typedef struct Query{
    int source;
    int target;
    long long binary_distance;
    long long ternary_distance;
    long long fib_distance;
    long long binomial_distance;
}Query;

#define QUERY_COUNT 1000

static void free_graph(Graph *g){
    free(g->offset);
    free(g->edges);
    g->n=0;
    g->m=0;
    g->offset=NULL;
    g->edges=NULL;
}

static int read_header(const char *path,Graph *g){
    FILE *fp=fopen(path,"r");
    if(fp==NULL){
        perror(path);
        return 0;
    }

    char line[256];
    char format[16];
    int found=0;

    while(fgets(line,sizeof(line),fp)!=NULL){
        if(sscanf(line," p %15s %d %d",format,&g->n,&g->m)==3){
            found=1;
            break;
        }
    }

    fclose(fp);

    if(!found||strcmp(format,"sp")!=0||g->n<=0||g->m<0){
        fprintf(stderr,"Invalid DIMACS graph header\n");
        return 0;
    }

    return 1;
}

static int count_degree(const char *path,const Graph *g,int *degree){
    FILE *fp=fopen(path,"r");
    if(fp==NULL){
        perror(path);
        return 0;
    }

    char line[256];
    char tag;
    int u,v,w;
    int count=0;
    int ok=1;

    while(fgets(line,sizeof(line),fp)!=NULL){
        if(sscanf(line," %c",&tag)!=1||tag!='a'){
            continue;
        }

        if(sscanf(line," a %d %d %d",&u,&v,&w)!=3||
           u<1||u>g->n||v<1||v>g->n||w<0){
            fprintf(stderr,"Invalid arc or negative edge weight\n");
            ok=0;
            break;
        }

        degree[u-1]++;
        count++;
    }

    fclose(fp);

    if(ok&&count!=g->m){
        fprintf(stderr,"Arc count mismatch: expected %d, read %d\n",
                g->m,count);
        ok=0;
    }

    return ok;
}

static int fill_edges(const char *path,Graph *g,int *cursor){
    FILE *fp=fopen(path,"r");
    if(fp==NULL){
        perror(path);
        return 0;
    }

    char line[256];
    char tag;
    int u,v,w;
    int count=0;
    int ok=1;

    while(fgets(line,sizeof(line),fp)!=NULL){
        if(sscanf(line," %c",&tag)!=1||tag!='a'){
            continue;
        }

        if(sscanf(line," a %d %d %d",&u,&v,&w)!=3||
           u<1||u>g->n||v<1||v>g->n||w<0){
            fprintf(stderr,"Invalid arc or negative edge weight\n");
            ok=0;
            break;
        }

        int source=u-1;
        int pos=cursor[source];

        if(pos>=g->offset[source+1]){
            fprintf(stderr,"Out-degree changed between file scans\n");
            ok=0;
            break;
        }

        g->edges[pos].t=v-1;
        g->edges[pos].w=w;
        cursor[source]++;
        count++;
    }

    fclose(fp);

    if(ok&&count!=g->m){
        fprintf(stderr,"Arc count mismatch during second scan\n");
        ok=0;
    }

    if(ok){
        for(int i=0;i<g->n;i++){
            if(cursor[i]!=g->offset[i+1]){
                fprintf(stderr,"CSR edge range was not filled completely\n");
                ok=0;
                break;
            }
        }
    }

    return ok;
}

static int load_graph(const char *path,Graph *g){
    int *degree=NULL;
    int *cursor=NULL;
    int ok=0;

    g->n=0;
    g->m=0;
    g->offset=NULL;
    g->edges=NULL;

    if(!read_header(path,g)){
        goto cleanup;
    }

    degree=(int*)calloc((size_t)g->n,sizeof(int));
    if(degree==NULL){
        fprintf(stderr,"Cannot allocate degree array\n");
        goto cleanup;
    }

    if(!count_degree(path,g,degree)){
        goto cleanup;
    }

    g->offset=(int*)malloc(((size_t)g->n+1)*sizeof(int));
    if(g->offset==NULL){
        fprintf(stderr,"Cannot allocate offset array\n");
        goto cleanup;
    }

    g->offset[0]=0;
    for(int i=0;i<g->n;i++){
        if(g->offset[i]>INT_MAX-degree[i]){
            fprintf(stderr,"Arc count exceeds int range\n");
            goto cleanup;
        }
        g->offset[i+1]=g->offset[i]+degree[i];
    }

    if(g->offset[g->n]!=g->m){
        fprintf(stderr,"Total arc count does not match the header\n");
        goto cleanup;
    }

    if(g->m>0){
        g->edges=(road*)malloc((size_t)g->m*sizeof(road));
        if(g->edges==NULL){
            fprintf(stderr,"Cannot allocate edge array\n");
            goto cleanup;
        }
    }

    cursor=(int*)malloc((size_t)g->n*sizeof(int));
    if(cursor==NULL){
        fprintf(stderr,"Cannot allocate cursor array\n");
        goto cleanup;
    }

    memcpy(cursor,g->offset,(size_t)g->n*sizeof(int));

    if(!fill_edges(path,g,cursor)){
        goto cleanup;
    }

    ok=1;

cleanup:
    free(degree);
    free(cursor);

    if(!ok){
        free_graph(g);
    }

    return ok;
}

static int heap_item_less(const HeapItem *a,const HeapItem *b){
    METRIC_ADD(comparisons,1);
    if(a->key!=b->key){
        return a->key<b->key;
    }
    return a->vertex<b->vertex;
}

static void binary_swap(BinaryHeap *heap,int i,int j){
    METRIC_ADD(array_swaps,1);
    HeapItem temp=heap->items[i];
    heap->items[i]=heap->items[j];
    heap->items[j]=temp;

    heap->position[heap->items[i].vertex]=i;
    heap->position[heap->items[j].vertex]=j;
}

static void binary_sift_up(BinaryHeap *heap,int i){
    while(i>0){
        METRIC_ADD(upward_steps,1);
        int parent=(i-1)/2;
        if(!heap_item_less(&heap->items[i],&heap->items[parent])){
            break;
        }
        binary_swap(heap,i,parent);
        i=parent;
    }
}

static void binary_sift_down(BinaryHeap *heap,int i){
    while(1){
        METRIC_ADD(downward_iterations,1);
        int left=2*i+1;
        int right=left+1;
        int smallest=i;

        if(left<heap->size&&
           heap_item_less(&heap->items[left],&heap->items[smallest])){
            smallest=left;
        }
        if(right<heap->size&&
           heap_item_less(&heap->items[right],&heap->items[smallest])){
            smallest=right;
        }
        if(smallest==i){
            break;
        }

        binary_swap(heap,i,smallest);
        i=smallest;
    }
}

static int binary_heap_init(BinaryHeap *heap,int capacity){
    heap->items=NULL;
    heap->position=NULL;
    heap->size=0;
    heap->capacity=capacity;

    if(capacity<=0){
        return 0;
    }

    heap->items=(HeapItem*)malloc((size_t)capacity*sizeof(HeapItem));
    heap->position=(int*)malloc((size_t)capacity*sizeof(int));

    if(heap->items==NULL||heap->position==NULL){
        free(heap->items);
        free(heap->position);
        heap->items=NULL;
        heap->position=NULL;
        return 0;
    }

    for(int i=0;i<capacity;i++){
        heap->position[i]=-1;
    }

    return 1;
}

static void binary_heap_destroy(BinaryHeap *heap){
    free(heap->items);
    free(heap->position);
    heap->items=NULL;
    heap->position=NULL;
    heap->size=0;
    heap->capacity=0;
}

static int binary_heap_empty(void *data){
    BinaryHeap *heap=(BinaryHeap*)data;
    return heap->size==0;
}

static int binary_heap_insert(void *data,int vertex,long long key){
    BinaryHeap *heap=(BinaryHeap*)data;

    if(vertex<0||vertex>=heap->capacity||
       heap->size>=heap->capacity||heap->position[vertex]!=-1){
        return 0;
    }

    int i=heap->size++;
    METRIC_MAX(max_queue_size,heap->size);
    heap->items[i].vertex=vertex;
    heap->items[i].key=key;
    heap->position[vertex]=i;
    binary_sift_up(heap,i);
    return 1;
}

static int binary_heap_extract_min(void *data,HeapItem *result){
    BinaryHeap *heap=(BinaryHeap*)data;
    if(heap->size==0){
        return 0;
    }

    *result=heap->items[0];
    heap->position[result->vertex]=-1;
    heap->size--;

    if(heap->size>0){
        heap->items[0]=heap->items[heap->size];
        heap->position[heap->items[0].vertex]=0;
        binary_sift_down(heap,0);
    }

    return 1;
}

static int binary_heap_decrease_key(void *data,int vertex,long long new_key){
    BinaryHeap *heap=(BinaryHeap*)data;

    if(vertex<0||vertex>=heap->capacity){
        return 0;
    }

    int i=heap->position[vertex];
    if(i<0||new_key>heap->items[i].key){
        return 0;
    }

    heap->items[i].key=new_key;
    binary_sift_up(heap,i);
    return 1;
}



static void ternary_swap(TernaryHeap *heap,int i,int j){
    METRIC_ADD(array_swaps,1);
    HeapItem temp=heap->items[i];
    heap->items[i]=heap->items[j];
    heap->items[j]=temp;

    heap->position[heap->items[i].vertex]=i;
    heap->position[heap->items[j].vertex]=j;
}

/* In a 0-based 3-ary heap, parent(i) is (i-1)/3. */
static void ternary_sift_up(TernaryHeap *heap,int i){
    while(i>0){
        METRIC_ADD(upward_steps,1);
        int parent=(i-1)/3;
        if(!heap_item_less(&heap->items[i],&heap->items[parent])){
            break;
        }
        ternary_swap(heap,i,parent);
        i=parent;
    }
}

/* Children of i are 3*i+1, 3*i+2, and 3*i+3. The heap never holds more
   than n entries, and 3*n stays far inside int range for these graphs,
   so 32-bit child arithmetic is safe and avoids 64-bit index math. */
static void ternary_sift_down(TernaryHeap *heap,int i){
    while(1){
        METRIC_ADD(downward_iterations,1);
        int first_child=3*i+1;
        int smallest=i;

        if(first_child<heap->size&&
           heap_item_less(&heap->items[first_child],
                          &heap->items[smallest])){
            smallest=first_child;
        }
        if(first_child+1<heap->size&&
           heap_item_less(&heap->items[first_child+1],
                          &heap->items[smallest])){
            smallest=first_child+1;
        }
        if(first_child+2<heap->size&&
           heap_item_less(&heap->items[first_child+2],
                          &heap->items[smallest])){
            smallest=first_child+2;
        }
        if(smallest==i){
            break;
        }

        ternary_swap(heap,i,smallest);
        i=smallest;
    }
}

static int ternary_heap_init(TernaryHeap *heap,int capacity){
    heap->items=NULL;
    heap->position=NULL;
    heap->size=0;
    heap->capacity=capacity;

    if(capacity<=0){
        return 0;
    }

    heap->items=(HeapItem*)malloc((size_t)capacity*sizeof(HeapItem));
    heap->position=(int*)malloc((size_t)capacity*sizeof(int));

    if(heap->items==NULL||heap->position==NULL){
        free(heap->items);
        free(heap->position);
        heap->items=NULL;
        heap->position=NULL;
        return 0;
    }

    for(int i=0;i<capacity;i++){
        heap->position[i]=-1;
    }

    return 1;
}

static void ternary_heap_destroy(TernaryHeap *heap){
    free(heap->items);
    free(heap->position);
    heap->items=NULL;
    heap->position=NULL;
    heap->size=0;
    heap->capacity=0;
}

static int ternary_heap_empty(void *data){
    TernaryHeap *heap=(TernaryHeap*)data;
    return heap->size==0;
}

static int ternary_heap_insert(void *data,int vertex,long long key){
    TernaryHeap *heap=(TernaryHeap*)data;

    if(vertex<0||vertex>=heap->capacity||
       heap->size>=heap->capacity||heap->position[vertex]!=-1){
        return 0;
    }

    int i=heap->size++;
    METRIC_MAX(max_queue_size,heap->size);
    heap->items[i].vertex=vertex;
    heap->items[i].key=key;
    heap->position[vertex]=i;
    ternary_sift_up(heap,i);
    return 1;
}

static int ternary_heap_extract_min(void *data,HeapItem *result){
    TernaryHeap *heap=(TernaryHeap*)data;
    if(heap->size==0){
        return 0;
    }

    *result=heap->items[0];
    heap->position[result->vertex]=-1;
    heap->size--;

    if(heap->size>0){
        heap->items[0]=heap->items[heap->size];
        heap->position[heap->items[0].vertex]=0;
        ternary_sift_down(heap,0);
    }

    return 1;
}

static int ternary_heap_decrease_key(void *data,int vertex,
                                     long long new_key){
    TernaryHeap *heap=(TernaryHeap*)data;

    if(vertex<0||vertex>=heap->capacity){
        return 0;
    }

    int i=heap->position[vertex];
    if(i<0||new_key>heap->items[i].key){
        return 0;
    }

    heap->items[i].key=new_key;
    ternary_sift_up(heap,i);
    return 1;
}

static int binomial_node_less(const BinomialNode *a,
                              const BinomialNode *b){
    METRIC_ADD(comparisons,1);
    if(a->key!=b->key){
        return a->key<b->key;
    }
    return a->vertex<b->vertex;
}

static BinomialNode *binomial_merge_roots(BinomialNode *a,
                                          BinomialNode *b){
    BinomialNode dummy;
    BinomialNode *tail=&dummy;
    dummy.sibling=NULL;
#ifdef COUNT_HEAP_OPS
    int merged_root_count=0;
#endif

    while(a!=NULL&&b!=NULL){
        METRIC_ADD(root_list_steps,1);
        if(a->degree<=b->degree){
            tail->sibling=a;
            a=a->sibling;
        }else{
            tail->sibling=b;
            b=b->sibling;
        }
        tail=tail->sibling;
#ifdef COUNT_HEAP_OPS
        merged_root_count++;
#endif
    }

    tail->sibling=(a!=NULL)?a:b;
#ifdef COUNT_HEAP_OPS
    /* Counting-only pass over the merged list. Keeping it behind the
       counter switch removes a full root-list traversal from every
       insert and every extract-min in the timing build. */
    {
        BinomialNode *rest=tail->sibling;
        while(rest!=NULL){
            merged_root_count++;
            rest=rest->sibling;
        }
        METRIC_MAX(max_root_count,merged_root_count);
    }
#endif
    return dummy.sibling;
}

static void binomial_link(BinomialNode *child,BinomialNode *parent){
    METRIC_ADD(tree_links,1);
    child->parent=parent;
    child->sibling=parent->child;
    parent->child=child;
    parent->degree++;
    METRIC_MAX(max_degree,parent->degree);
}

static void binomial_union_roots(BinomialHeap *heap,
                                 BinomialNode *other_roots){
    BinomialNode *head=binomial_merge_roots(heap->roots,other_roots);
    if(head==NULL){
        heap->roots=NULL;
        return;
    }

    BinomialNode *previous=NULL;
    BinomialNode *current=head;
    BinomialNode *next=current->sibling;

    while(next!=NULL){
        METRIC_ADD(root_list_steps,1);
        if(current->degree!=next->degree||
           (next->sibling!=NULL&&
            next->sibling->degree==current->degree)){
            previous=current;
            current=next;
        }else if(binomial_node_less(current,next)){
            current->sibling=next->sibling;
            binomial_link(next,current);
        }else{
            if(previous==NULL){
                head=next;
            }else{
                previous->sibling=next;
            }
            binomial_link(current,next);
            current=next;
        }

        next=current->sibling;
    }

    heap->roots=head;
}

static int binomial_heap_init(BinomialHeap *heap,int n){
    heap->nodes=NULL;
    heap->node_of=NULL;
    heap->roots=NULL;
    heap->n=n;
    heap->size=0;
    heap->node_count=0;

    if(n<=0){
        return 0;
    }

    heap->nodes=(BinomialNode*)calloc((size_t)n,sizeof(BinomialNode));
    heap->node_of=(BinomialNode**)calloc((size_t)n,sizeof(BinomialNode*));

    if(heap->nodes==NULL||heap->node_of==NULL){
        free(heap->nodes);
        free(heap->node_of);
        heap->nodes=NULL;
        heap->node_of=NULL;
        return 0;
    }

    return 1;
}

static void binomial_heap_destroy(BinomialHeap *heap){
    free(heap->nodes);
    free(heap->node_of);
    heap->nodes=NULL;
    heap->node_of=NULL;
    heap->roots=NULL;
    heap->n=0;
    heap->size=0;
    heap->node_count=0;
}

static int binomial_heap_empty(void *data){
    BinomialHeap *heap=(BinomialHeap*)data;
    return heap->size==0;
}

static int binomial_heap_insert(void *data,int vertex,long long key){
    BinomialHeap *heap=(BinomialHeap*)data;

    if(vertex<0||vertex>=heap->n||heap->size>=heap->n||
       heap->node_count>=heap->n||heap->node_of[vertex]!=NULL){
        return 0;
    }

    BinomialNode *node=&heap->nodes[heap->node_count++];
    node->vertex=vertex;
    node->key=key;
    node->degree=0;
    node->parent=NULL;
    node->child=NULL;
    node->sibling=NULL;
    heap->node_of[vertex]=node;

    binomial_union_roots(heap,node);
    heap->size++;
    METRIC_MAX(max_queue_size,heap->size);
    return 1;
}

static int binomial_heap_extract_min(void *data,HeapItem *result){
    BinomialHeap *heap=(BinomialHeap*)data;
    if(heap->roots==NULL||heap->size==0){
        return 0;
    }

    BinomialNode *min_node=heap->roots;
    BinomialNode *min_previous=NULL;
    BinomialNode *previous=heap->roots;
    BinomialNode *current=heap->roots->sibling;

    while(current!=NULL){
        METRIC_ADD(root_list_steps,1);
        if(binomial_node_less(current,min_node)){
            min_node=current;
            min_previous=previous;
        }
        previous=current;
        current=current->sibling;
    }

    if(min_previous==NULL){
        heap->roots=min_node->sibling;
    }else{
        min_previous->sibling=min_node->sibling;
    }

    BinomialNode *child=min_node->child;
    BinomialNode *reversed_children=NULL;

    while(child!=NULL){
        METRIC_ADD(children_promoted,1);
        BinomialNode *next=child->sibling;
        child->parent=NULL;
        child->sibling=reversed_children;
        reversed_children=child;
        child=next;
    }

    result->vertex=min_node->vertex;
    result->key=min_node->key;
    heap->node_of[result->vertex]=NULL;
    heap->size--;

    min_node->parent=NULL;
    min_node->child=NULL;
    min_node->sibling=NULL;
    min_node->degree=0;

    binomial_union_roots(heap,reversed_children);
    return 1;
}

static int binomial_heap_decrease_key(void *data,int vertex,
                                      long long new_key){
    BinomialHeap *heap=(BinomialHeap*)data;

    if(vertex<0||vertex>=heap->n||heap->node_of[vertex]==NULL){
        return 0;
    }

    BinomialNode *node=heap->node_of[vertex];
    if(new_key>node->key){
        return 0;
    }

    node->key=new_key;

    while(node->parent!=NULL){
        METRIC_ADD(upward_steps,1);
        if(!binomial_node_less(node,node->parent)){
            break;
        }
        BinomialNode *parent=node->parent;
        int vertex_temp=node->vertex;
        long long key_temp=node->key;

        METRIC_ADD(key_swaps,1);
        node->vertex=parent->vertex;
        node->key=parent->key;
        parent->vertex=vertex_temp;
        parent->key=key_temp;

        heap->node_of[node->vertex]=node;
        heap->node_of[parent->vertex]=parent;
        node=parent;
    }

    return 1;
}

static int fib_node_less(const FibNode *a,const FibNode *b){
    METRIC_ADD(comparisons,1);
    if(a->key!=b->key){
        return a->key<b->key;
    }
    return a->vertex<b->vertex;
}

static void fib_remove_from_list(FibNode *x){
    x->left->right=x->right;
    x->right->left=x->left;
    x->left=x;
    x->right=x;
}

static void fib_add_root(FibHeap *heap,FibNode *x){
    x->parent=NULL;
    x->mark=0;

    if(heap->min==NULL){
        x->left=x;
        x->right=x;
        heap->min=x;
        METRIC_MAX(max_root_count,1);
        return;
    }

    x->left=heap->min;
    x->right=heap->min->right;
    heap->min->right->left=x;
    heap->min->right=x;

    if(fib_node_less(x,heap->min)){
        heap->min=x;
    }
}

static void fib_link(FibNode *child,FibNode *parent){
    METRIC_ADD(tree_links,1);
    fib_remove_from_list(child);
    child->parent=parent;
    child->mark=0;

    if(parent->child==NULL){
        child->left=child;
        child->right=child;
        parent->child=child;
    }else{
        child->left=parent->child;
        child->right=parent->child->right;
        parent->child->right->left=child;
        parent->child->right=child;
    }

    parent->degree++;
    METRIC_MAX(max_degree,parent->degree);
}

static int fib_consolidate(FibHeap *heap){
    if(heap->min==NULL){
        return 1;
    }

    int root_count=0;
    FibNode *current=heap->min;

    do{
        if(root_count>=heap->n){
            fprintf(stderr,"Invalid Fibonacci heap root list\n");
            return 0;
        }
        heap->root_buffer[root_count++]=current;
        METRIC_ADD(root_list_steps,1);
        current=current->right;
    }while(current!=heap->min);
    METRIC_MAX(max_root_count,root_count);

    /* The degree table stays all-NULL between consolidations: entries
       consumed during linking are cleared immediately, and surviving
       entries are cleared below as they move back to the root list.
       Only the touched range [0,hi] is scanned, which removes the former
       O(degree_capacity) memset and full-table scan per extract-min. */
    int hi=-1;

    for(int i=0;i<root_count;i++){
        FibNode *x=heap->root_buffer[i];

        if(x->parent!=NULL){
            continue;
        }

        int degree=x->degree;
        while(degree<heap->degree_capacity&&
              heap->degree_table[degree]!=NULL){
            FibNode *y=heap->degree_table[degree];

            if(fib_node_less(y,x)){
                FibNode *temp=x;
                x=y;
                y=temp;
            }

            heap->degree_table[degree]=NULL;
            fib_link(y,x);
            degree=x->degree;
        }

        if(degree>=heap->degree_capacity){
            fprintf(stderr,"Fibonacci heap degree bound exceeded\n");
            return 0;
        }

        heap->degree_table[degree]=x;
        if(degree>hi){
            hi=degree;
        }
    }

    heap->min=NULL;

    for(int i=0;i<=hi;i++){
        FibNode *x=heap->degree_table[i];
        if(x==NULL){
            continue;
        }

        heap->degree_table[i]=NULL;
        x->left=x;
        x->right=x;
        x->parent=NULL;
        x->mark=0;
        fib_add_root(heap,x);
    }

    return 1;
}

static int fib_heap_init(FibHeap *heap,int n){
    heap->nodes=NULL;
    heap->root_buffer=NULL;
    heap->degree_table=NULL;
    heap->degree_capacity=0;
    heap->n=n;
    heap->size=0;
    heap->min=NULL;

    if(n<=0){
        return 0;
    }

    int capacity=3;
    for(int x=n;x>1;x>>=1){
        capacity+=2;
    }
    heap->degree_capacity=capacity;

    heap->nodes=(FibNode*)calloc((size_t)n,sizeof(FibNode));
    heap->root_buffer=(FibNode**)malloc((size_t)n*sizeof(FibNode*));
    heap->degree_table=(FibNode**)calloc((size_t)capacity,sizeof(FibNode*));

    if(heap->nodes==NULL||heap->root_buffer==NULL||
       heap->degree_table==NULL){
        free(heap->nodes);
        free(heap->root_buffer);
        free(heap->degree_table);
        heap->nodes=NULL;
        heap->root_buffer=NULL;
        heap->degree_table=NULL;
        heap->degree_capacity=0;
        return 0;
    }

    for(int i=0;i<n;i++){
        heap->nodes[i].vertex=i;
    }

    return 1;
}

static void fib_heap_destroy(FibHeap *heap){
    free(heap->nodes);
    free(heap->root_buffer);
    free(heap->degree_table);
    heap->nodes=NULL;
    heap->root_buffer=NULL;
    heap->degree_table=NULL;
    heap->degree_capacity=0;
    heap->n=0;
    heap->size=0;
    heap->min=NULL;
}

static int fib_heap_empty(void *data){
    FibHeap *heap=(FibHeap*)data;
    return heap->size==0;
}

static int fib_heap_insert(void *data,int vertex,long long key){
    FibHeap *heap=(FibHeap*)data;

    if(vertex<0||vertex>=heap->n||
       heap->size>=heap->n||heap->nodes[vertex].in_heap){
        return 0;
    }

    FibNode *x=&heap->nodes[vertex];
    x->key=key;
    x->degree=0;
    x->mark=0;
    x->in_heap=1;
    x->parent=NULL;
    x->child=NULL;
    x->left=x;
    x->right=x;

    fib_add_root(heap,x);
    heap->size++;
    METRIC_MAX(max_queue_size,heap->size);
    return 1;
}

static void fib_cut(FibHeap *heap,FibNode *x,FibNode *parent){
    METRIC_ADD(cuts,1);
    if(parent->child==x){
        if(x->right==x){
            parent->child=NULL;
        }else{
            parent->child=x->right;
        }
    }

    fib_remove_from_list(x);
    parent->degree--;
    fib_add_root(heap,x);
}

static void fib_cascading_cut(FibHeap *heap,FibNode *x){
    while(x!=NULL){
        METRIC_ADD(cascade_steps,1);
        FibNode *parent=x->parent;
        if(parent==NULL){
            x->mark=0;
            return;
        }

        if(!x->mark){
            x->mark=1;
            return;
        }

        fib_cut(heap,x,parent);
        x=parent;
    }
}

static int fib_heap_decrease_key(void *data,int vertex,long long new_key){
    FibHeap *heap=(FibHeap*)data;

    if(vertex<0||vertex>=heap->n){
        return 0;
    }

    FibNode *x=&heap->nodes[vertex];
    if(!x->in_heap||new_key>x->key){
        return 0;
    }

    x->key=new_key;
    FibNode *parent=x->parent;

    if(parent!=NULL&&fib_node_less(x,parent)){
        fib_cut(heap,x,parent);
        fib_cascading_cut(heap,parent);
    }

    if(heap->min==NULL||fib_node_less(x,heap->min)){
        heap->min=x;
    }

    return 1;
}

static int fib_heap_extract_min(void *data,HeapItem *result){
    FibHeap *heap=(FibHeap*)data;
    FibNode *z=heap->min;

    if(z==NULL){
        return 0;
    }

    result->vertex=z->vertex;
    result->key=z->key;

    while(z->child!=NULL){
        FibNode *child=z->child;
        METRIC_ADD(children_promoted,1);

        if(child->right==child){
            z->child=NULL;
        }else{
            z->child=child->right;
        }

        fib_remove_from_list(child);
        child->parent=NULL;
        child->mark=0;
        fib_add_root(heap,child);
    }

    z->degree=0;

    if(z->right==z){
        heap->min=NULL;
    }else{
        FibNode *next=z->right;
        fib_remove_from_list(z);
        heap->min=next;
    }

    heap->size--;
    z->in_heap=0;
    z->parent=NULL;
    z->child=NULL;
    z->left=z;
    z->right=z;
    z->mark=0;

    if(heap->min!=NULL&&!fib_consolidate(heap)){
        return 0;
    }

    return 1;
}

/* Per-query cleanup. The Dijkstra loop below calls the concrete heap
   functions directly (via the generated specializations), so no vtable
   wrappers are needed for the hot operations. */
static int binary_pq_clear(void *data,const int *touched,int count){
    BinaryHeap *heap=(BinaryHeap*)data;

    for(int i=0;i<count;i++){
        int vertex=touched[i];
        if(vertex<0||vertex>=heap->capacity){
            return 0;
        }
        heap->position[vertex]=-1;
    }

    heap->size=0;
    return 1;
}

static int ternary_pq_clear(void *data,const int *touched,int count){
    TernaryHeap *heap=(TernaryHeap*)data;

    for(int i=0;i<count;i++){
        int vertex=touched[i];
        if(vertex<0||vertex>=heap->capacity){
            return 0;
        }
        heap->position[vertex]=-1;
    }

    heap->size=0;
    return 1;
}

static int fib_pq_clear(void *data,const int *touched,int count){
    FibHeap *heap=(FibHeap*)data;

    for(int i=0;i<count;i++){
        int vertex=touched[i];
        if(vertex<0||vertex>=heap->n){
            return 0;
        }
        heap->nodes[vertex].in_heap=0;
    }

    heap->min=NULL;
    heap->size=0;
    return 1;
}

static int binomial_pq_clear(void *data,const int *touched,int count){
    BinomialHeap *heap=(BinomialHeap*)data;

    for(int i=0;i<count;i++){
        int vertex=touched[i];
        if(vertex<0||vertex>=heap->n){
            return 0;
        }
        heap->node_of[vertex]=NULL;
    }

    heap->roots=NULL;
    heap->size=0;
    heap->node_count=0;
    return 1;
}


static int dijkstra_workspace_init(DijkstraWorkspace *workspace,int n){
    workspace->dist=NULL;
    workspace->state=NULL;
    workspace->touched=NULL;
    workspace->touched_count=0;
    workspace->n=n;

    if(n<=0){
        return 0;
    }

    workspace->dist=(long long*)malloc((size_t)n*sizeof(long long));
    workspace->state=(unsigned char*)calloc((size_t)n,sizeof(unsigned char));
    workspace->touched=(int*)malloc((size_t)n*sizeof(int));

    if(workspace->dist==NULL||workspace->state==NULL||
       workspace->touched==NULL){
        free(workspace->dist);
        free(workspace->state);
        free(workspace->touched);
        workspace->dist=NULL;
        workspace->state=NULL;
        workspace->touched=NULL;
        return 0;
    }

    for(int i=0;i<n;i++){
        workspace->dist[i]=DIST_INF;
    }

    return 1;
}

static void dijkstra_workspace_destroy(DijkstraWorkspace *workspace){
    free(workspace->dist);
    free(workspace->state);
    free(workspace->touched);
    workspace->dist=NULL;
    workspace->state=NULL;
    workspace->touched=NULL;
    workspace->touched_count=0;
    workspace->n=0;
}

/* One Dijkstra body compiled once per heap type. Passing the concrete
   heap functions as macro arguments lets the compiler inline insert,
   extract-min, decrease-key and the empty check. The previous version
   dispatched every operation through a function-pointer table, which
   blocked inlining on hundreds of millions of calls. */
#define DEFINE_DIJKSTRA_QUERY(NAME,HEAP_TYPE,EMPTY,INSERT,EXTRACT,DECREASE,CLEAR) \
static int dijkstra_query_##NAME(const Graph *g,int source,int target, \
                                 void *heap_void, \
                                 DijkstraWorkspace *workspace, \
                                 SearchStats *stats,long long *answer){ \
    HEAP_TYPE *heap=(HEAP_TYPE*)heap_void; \
    int ok=0; \
    \
    if(source<0||source>=g->n||target<0||target>=g->n|| \
       workspace->n!=g->n||stats==NULL|| \
       workspace->touched_count!=0||!EMPTY(heap)){ \
        fprintf(stderr,"Invalid point-to-point Dijkstra input\n"); \
        return 0; \
    } \
    \
    *answer=DIST_INF; \
    workspace->dist[source]=0; \
    workspace->state[source]=1; \
    workspace->touched[0]=source; \
    workspace->touched_count=1; \
    \
    if(!INSERT(heap,source,0)){ \
        fprintf(stderr,"Cannot insert source into priority queue\n"); \
        goto cleanup; \
    } \
    stats->pq_inserts++; \
    \
    while(!EMPTY(heap)){ \
        HeapItem item; \
        \
        if(!EXTRACT(heap,&item)){ \
            fprintf(stderr,"Priority queue extraction failed\n"); \
            goto cleanup; \
        } \
        stats->pq_extracts++; \
        \
        int u=item.vertex; \
        if(u<0||u>=g->n||workspace->state[u]!=1|| \
           item.key!=workspace->dist[u]){ \
            fprintf(stderr,"Priority queue returned an invalid item\n"); \
            goto cleanup; \
        } \
        \
        workspace->state[u]=2; \
        stats->settled_vertices++; \
        \
        if(u==target){ \
            *answer=workspace->dist[u]; \
            ok=1; \
            goto cleanup; \
        } \
        \
        for(int i=g->offset[u];i<g->offset[u+1];i++){ \
            stats->scanned_arcs++; \
            int v=g->edges[i].t; \
            int weight=g->edges[i].w; \
            \
            if(v<0||v>=g->n||weight<0){ \
                fprintf(stderr,"Invalid edge in graph\n"); \
                goto cleanup; \
            } \
            \
            if(workspace->dist[u]>DIST_INF-(long long)weight){ \
                continue; \
            } \
            \
            long long candidate=workspace->dist[u]+(long long)weight; \
            \
            if(candidate<workspace->dist[v]){ \
                if(workspace->state[v]==2){ \
                    fprintf(stderr, \
                            "A finalized distance was decreased\n"); \
                    goto cleanup; \
                } \
                \
                if(workspace->state[v]==0){ \
                    if(workspace->touched_count>=workspace->n){ \
                        fprintf(stderr,"Too many touched vertices\n"); \
                        goto cleanup; \
                    } \
                    \
                    workspace->touched[workspace->touched_count++]=v; \
                    workspace->dist[v]=candidate; \
                    workspace->state[v]=1; \
                    \
                    if(!INSERT(heap,v,candidate)){ \
                        fprintf(stderr, \
                                "Priority queue insertion failed\n"); \
                        goto cleanup; \
                    } \
                    stats->pq_inserts++; \
                }else{ \
                    workspace->dist[v]=candidate; \
                    if(!DECREASE(heap,v,candidate)){ \
                        fprintf(stderr, \
                                "Priority queue decrease-key failed\n"); \
                        goto cleanup; \
                    } \
                    stats->pq_decrease_keys++; \
                } \
            } \
        } \
    } \
    \
    ok=1; \
    \
cleanup: \
    if(!CLEAR(heap,workspace->touched,workspace->touched_count)){ \
        fprintf(stderr,"Priority queue cleanup failed\n"); \
        ok=0; \
    } \
    \
    for(int i=0;i<workspace->touched_count;i++){ \
        int vertex=workspace->touched[i]; \
        workspace->dist[vertex]=DIST_INF; \
        workspace->state[vertex]=0; \
    } \
    workspace->touched_count=0; \
    \
    return ok; \
}

DEFINE_DIJKSTRA_QUERY(binary,BinaryHeap,
                      binary_heap_empty,binary_heap_insert,
                      binary_heap_extract_min,binary_heap_decrease_key,
                      binary_pq_clear)
DEFINE_DIJKSTRA_QUERY(ternary,TernaryHeap,
                      ternary_heap_empty,ternary_heap_insert,
                      ternary_heap_extract_min,ternary_heap_decrease_key,
                      ternary_pq_clear)
DEFINE_DIJKSTRA_QUERY(fibonacci,FibHeap,
                      fib_heap_empty,fib_heap_insert,
                      fib_heap_extract_min,fib_heap_decrease_key,
                      fib_pq_clear)
DEFINE_DIJKSTRA_QUERY(binomial,BinomialHeap,
                      binomial_heap_empty,binomial_heap_insert,
                      binomial_heap_extract_min,binomial_heap_decrease_key,
                      binomial_pq_clear)

typedef int (*QueryFn)(const Graph *,int,int,void *,
                       DijkstraWorkspace *,SearchStats *,long long *);

static unsigned int next_random(unsigned int *state){
    unsigned int x=*state;
    x^=x<<13;
    x^=x>>17;
    x^=x<<5;
    *state=x;
    return x;
}

static int generate_queries(Query *queries,int n,int count,
                            unsigned int seed){
    if(n<2){
        return 0;
    }
    /* A graph with n vertices has only n*(n-1) ordered pairs with
       source != target. Without this guard the deduplication loop below
       would spin forever on a small graph. */
    if((long long)n*(long long)(n-1)<(long long)count){
        fprintf(stderr,"Not enough distinct source-target pairs\n");
        return 0;
    }
    if(seed==0){
        seed=1;
    }

    for(int i=0;i<count;i++){
        int duplicate;
        do{
            int source=(int)(next_random(&seed)%(unsigned int)n);
            int target=(int)(next_random(&seed)%(unsigned int)n);

            if(target==source){
                target=(target+1)%n;
            }

            duplicate=0;
            for(int j=0;j<i;j++){
                if(queries[j].source==source&&
                   queries[j].target==target){
                    duplicate=1;
                    break;
                }
            }

            if(!duplicate){
                queries[i].source=source;
                queries[i].target=target;
                queries[i].binary_distance=DIST_INF;
                queries[i].ternary_distance=DIST_INF;
                queries[i].fib_distance=DIST_INF;
                queries[i].binomial_distance=DIST_INF;
            }
        }while(duplicate);
    }

    return 1;
}

static int measure_query_batch(const Graph *g,Query *queries,int count,
                               QueryFn query_fn,void *heap,
                               DijkstraWorkspace *workspace,
                               int heap_type,double *seconds,
                               SearchStats *stats){
    memset(stats,0,sizeof(*stats));
    clock_t start=clock();
    if(start==(clock_t)-1){
        fprintf(stderr,"clock() failed\n");
        return 0;
    }

    for(int i=0;i<count;i++){
        long long answer;
        if(!query_fn(g,queries[i].source,queries[i].target,
                     heap,workspace,stats,&answer)){
            return 0;
        }

        if(heap_type==0){
            queries[i].binary_distance=answer;
        }else if(heap_type==1){
            queries[i].fib_distance=answer;
        }else if(heap_type==2){
            queries[i].binomial_distance=answer;
        }else{
            queries[i].ternary_distance=answer;
        }
    }

    clock_t finish=clock();
    if(finish==(clock_t)-1){
        fprintf(stderr,"clock() failed\n");
        return 0;
    }

    *seconds=(double)(finish-start)/(double)CLOCKS_PER_SEC;
    return 1;
}

static int run_binary_batch(const Graph *g,Query *queries,int count,
                            double *seconds,SearchStats *stats){
    BinaryHeap heap;
    DijkstraWorkspace workspace;

    if(!binary_heap_init(&heap,g->n)){
        fprintf(stderr,"Cannot initialize binary heap\n");
        return 0;
    }

    if(!dijkstra_workspace_init(&workspace,g->n)){
        fprintf(stderr,"Cannot initialize binary heap workspace\n");
        binary_heap_destroy(&heap);
        return 0;
    }

    int ok=measure_query_batch(g,queries,count,dijkstra_query_binary,
                               &heap,&workspace,0,seconds,stats);

    dijkstra_workspace_destroy(&workspace);
    binary_heap_destroy(&heap);
    return ok;
}

static int run_fib_batch(const Graph *g,Query *queries,int count,
                         double *seconds,SearchStats *stats){
    FibHeap heap;
    DijkstraWorkspace workspace;

    if(!fib_heap_init(&heap,g->n)){
        fprintf(stderr,"Cannot initialize Fibonacci heap\n");
        return 0;
    }

    if(!dijkstra_workspace_init(&workspace,g->n)){
        fprintf(stderr,"Cannot initialize Fibonacci heap workspace\n");
        fib_heap_destroy(&heap);
        return 0;
    }

    int ok=measure_query_batch(g,queries,count,dijkstra_query_fibonacci,
                               &heap,&workspace,1,seconds,stats);

    dijkstra_workspace_destroy(&workspace);
    fib_heap_destroy(&heap);
    return ok;
}


static int run_binomial_batch(const Graph *g,Query *queries,int count,
                              double *seconds,SearchStats *stats){
    BinomialHeap heap;

    if(!binomial_heap_init(&heap,g->n)){
        fprintf(stderr,"Cannot initialize binomial heap\n");
        return 0;
    }

    DijkstraWorkspace workspace;
    if(!dijkstra_workspace_init(&workspace,g->n)){
        fprintf(stderr,"Cannot initialize binomial heap workspace\n");
        binomial_heap_destroy(&heap);
        return 0;
    }

    int ok=measure_query_batch(g,queries,count,dijkstra_query_binomial,
                               &heap,&workspace,2,seconds,stats);

    dijkstra_workspace_destroy(&workspace);
    binomial_heap_destroy(&heap);
    return ok;
}

static int run_ternary_batch(const Graph *g,Query *queries,int count,
                             double *seconds,SearchStats *stats){
    TernaryHeap heap;
    DijkstraWorkspace workspace;

    if(!ternary_heap_init(&heap,g->n)){
        fprintf(stderr,"Cannot initialize ternary heap\n");
        return 0;
    }

    if(!dijkstra_workspace_init(&workspace,g->n)){
        fprintf(stderr,"Cannot initialize ternary heap workspace\n");
        ternary_heap_destroy(&heap);
        return 0;
    }

    int ok=measure_query_batch(g,queries,count,dijkstra_query_ternary,
                               &heap,&workspace,3,seconds,stats);

    dijkstra_workspace_destroy(&workspace);
    ternary_heap_destroy(&heap);
    return ok;
}

static long long csv_distance(long long distance){
    return distance==DIST_INF?-1:distance;
}


typedef int (*BatchRunner)(const Graph *,Query *,int,double *,SearchStats *);

static int benchmark_dataset(const char *name,const char *path,int index,
                             FILE *result_file,FILE *query_file,
                             FILE *operation_file){
    Graph g={0,0,NULL,NULL};
    Query *queries=NULL;
    int ok=0;
    unsigned int seed=0;
    const char *heap_names[4]={
        "binary","ternary","fibonacci","binomial"
    };
    BatchRunner runners[4]={
        run_binary_batch,run_ternary_batch,run_fib_batch,
        run_binomial_batch
    };
    double seconds[4]={0.0,0.0,0.0,0.0};
    SearchStats stats[4]={{0}};
    int matching=0;

    if(!load_graph(path,&g)){
        goto cleanup;
    }

    if(g.n<2){
        fprintf(stderr,"Graph %s has too few vertices\n",name);
        goto cleanup;
    }

    queries=(Query*)malloc((size_t)QUERY_COUNT*sizeof(Query));
    if(queries==NULL){
        fprintf(stderr,"Cannot allocate query array\n");
        goto cleanup;
    }

    seed=0x9e3779b9u^(unsigned int)g.n;
    if(!generate_queries(queries,g.n,QUERY_COUNT,seed)){
        fprintf(stderr,"Cannot generate query pairs\n");
        goto cleanup;
    }

    for(int order=0;order<4;order++){
        int heap_type=(order+index)%4;
        if(!runners[heap_type](&g,queries,QUERY_COUNT,
                               &seconds[heap_type],&stats[heap_type])){
            goto cleanup;
        }
    }

    for(int i=0;i<QUERY_COUNT;i++){
        long long answer=queries[i].binary_distance;
        if(answer==queries[i].ternary_distance&&
           answer==queries[i].fib_distance&&
           answer==queries[i].binomial_distance){
            matching++;
        }
    }

    for(int i=0;i<QUERY_COUNT;i++){
        if(fprintf(query_file,"%s,%d,%d,%d,%lld,%lld,%lld,%lld\n",
                   name,i+1,queries[i].source+1,queries[i].target+1,
                   csv_distance(queries[i].binary_distance),
                   csv_distance(queries[i].ternary_distance),
                   csv_distance(queries[i].fib_distance),
                   csv_distance(queries[i].binomial_distance))<0){
            fprintf(stderr,"Cannot write query CSV\n");
            goto cleanup;
        }
    }

    if(fprintf(result_file,
               "%s,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,"
               "%.3f,%.3f,%.3f,%.3f,%d\n",
               name,g.n,g.m,QUERY_COUNT,seconds[0],seconds[1],
               seconds[2],seconds[3],
               seconds[0]*1000.0/QUERY_COUNT,
               seconds[1]*1000.0/QUERY_COUNT,
               seconds[2]*1000.0/QUERY_COUNT,
               seconds[3]*1000.0/QUERY_COUNT,matching)<0){
        fprintf(stderr,"Cannot write result CSV\n");
        goto cleanup;
    }

    for(int h=0;h<4;h++){
        if(fprintf(operation_file,
                   "%s,%s,%d,%lld,%lld,%lld,%lld,%lld\n",
                   name,heap_names[h],QUERY_COUNT,
                   stats[h].scanned_arcs,stats[h].settled_vertices,
                   stats[h].pq_inserts,stats[h].pq_extracts,
                   stats[h].pq_decrease_keys)<0){
            fprintf(stderr,"Cannot write operation CSV\n");
            goto cleanup;
        }
    }

    printf("%s: %d vertices, %d arcs, %d queries\n",
           name,g.n,g.m,QUERY_COUNT);
    printf("  Binary heap:    %.3f s (%.3f ms/query)\n",
           seconds[0],seconds[0]*1000.0/QUERY_COUNT);
    printf("  Ternary heap:   %.3f s (%.3f ms/query)\n",
           seconds[1],seconds[1]*1000.0/QUERY_COUNT);
    printf("  Fibonacci heap: %.3f s (%.3f ms/query)\n",
           seconds[2],seconds[2]*1000.0/QUERY_COUNT);
    printf("  Binomial heap:  %.3f s (%.3f ms/query)\n",
           seconds[3],seconds[3]*1000.0/QUERY_COUNT);
    printf("  Matching answers: %d/%d\n",matching,QUERY_COUNT);
    fflush(stdout);

    ok=(matching==QUERY_COUNT);

cleanup:
    free(queries);
    free_graph(&g);
    return ok;
}


#ifdef COUNT_HEAP_OPS
static int collect_internal_stats_dataset(const char *name,const char *path,
                                          FILE *stats_file){
    Graph g={0,0,NULL,NULL};
    Query *queries=NULL;
    int ok=0;
    const char *heap_names[4]={
        "binary","ternary","fibonacci","binomial"
    };
    BatchRunner runners[4]={
        run_binary_batch,run_ternary_batch,run_fib_batch,
        run_binomial_batch
    };
    double ignored_seconds[4]={0.0,0.0,0.0,0.0};
    SearchStats work[4]={{0}};
    HeapInternalStats metrics[4]={{0}};
    int matching=0;
    unsigned int seed=0;

    if(!load_graph(path,&g)){
        goto cleanup;
    }

    queries=(Query*)malloc((size_t)QUERY_COUNT*sizeof(Query));
    if(queries==NULL){
        fprintf(stderr,"Cannot allocate query array\n");
        goto cleanup;
    }

    seed=0x9e3779b9u^(unsigned int)g.n;
    if(!generate_queries(queries,g.n,QUERY_COUNT,seed)){
        fprintf(stderr,"Cannot generate query pairs\n");
        goto cleanup;
    }

    for(int h=0;h<4;h++){
        printf("Counting %s / %s internal operations...\n",
               name,heap_names[h]);
        fflush(stdout);
#ifdef COUNT_HEAP_OPS
        active_heap_metrics=&metrics[h];
#endif
        int run_ok=runners[h](&g,queries,QUERY_COUNT,
                              &ignored_seconds[h],&work[h]);
#ifdef COUNT_HEAP_OPS
        active_heap_metrics=NULL;
#endif
        if(!run_ok){
            goto cleanup;
        }
    }

    for(int i=0;i<QUERY_COUNT;i++){
        long long answer=queries[i].binary_distance;
        if(answer==queries[i].ternary_distance&&
           answer==queries[i].fib_distance&&
           answer==queries[i].binomial_distance){
            matching++;
        }
    }

    for(int h=0;h<4;h++){
        if(fprintf(stats_file,
                   "%s,%s,%d,%lld,%lld,%lld,%lld,%lld,"
                   "%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,"
                   "%lld,%lld,%lld,%lld,%lld\n",
                   name,heap_names[h],QUERY_COUNT,
                   work[h].scanned_arcs,work[h].settled_vertices,
                   work[h].pq_inserts,work[h].pq_extracts,
                   work[h].pq_decrease_keys,
                   metrics[h].comparisons,metrics[h].upward_steps,
                   metrics[h].downward_iterations,
                   metrics[h].array_swaps,metrics[h].key_swaps,
                   metrics[h].tree_links,metrics[h].cuts,
                   metrics[h].cascade_steps,metrics[h].root_list_steps,
                   metrics[h].children_promoted,
                   metrics[h].max_queue_size,metrics[h].max_root_count,
                   metrics[h].max_degree)<0){
            fprintf(stderr,"Cannot write internal stats CSV\n");
            goto cleanup;
        }
    }

    printf("%s: matching answers %d/%d\n",name,matching,QUERY_COUNT);
    ok=(matching==QUERY_COUNT);

cleanup:
#ifdef COUNT_HEAP_OPS
    active_heap_metrics=NULL;
#endif
    free(queries);
    free_graph(&g);
    return ok;
}

static int run_internal_stats(void){
    const char *names[3]={"NY","COL","FLA"};
    const char *paths[3]={
        "data/dimacs/USA-road-d.NY.gr",
        "data/dimacs/USA-road-d.COL.gr",
        "data/dimacs/USA-road-d.FLA.gr"
    };

    FILE *stats_file=fopen("results/heap_internal_stats.csv","w");
    if(stats_file==NULL){
        perror("results/heap_internal_stats.csv");
        return 1;
    }

    fprintf(stats_file,
            "dataset,heap,queries,scanned_arcs,settled_vertices,"
            "pq_inserts,pq_extracts,pq_decrease_keys,comparisons,"
            "upward_steps,downward_iterations,array_swaps,key_swaps,"
            "tree_links,cuts,cascade_steps,root_list_steps,"
            "children_promoted,max_queue_size,max_root_count,max_degree\n");

    int ok=1;
    for(int i=0;i<3;i++){
        if(!collect_internal_stats_dataset(names[i],paths[i],stats_file)){
            ok=0;
            break;
        }
    }

    fclose(stats_file);
    if(!ok){
        fprintf(stderr,"Internal operation counting failed or answers differed\n");
        return 1;
    }

    printf("Saved results/heap_internal_stats.csv\n");
    return 0;
}
#endif

int main(int argc,char **argv){
    if(argc>1&&strcmp(argv[1],"--count-heap-ops")==0){
#ifdef COUNT_HEAP_OPS
        return run_internal_stats();
#else
        fprintf(stderr,
                "This binary was built without heap counters.\n"
                "Rebuild with -DCOUNT_HEAP_OPS to use --count-heap-ops.\n");
        return 1;
#endif
    }

    const char *names[3]={"NY","COL","FLA"};
    const char *paths[3]={
        "data/dimacs/USA-road-d.NY.gr",
        "data/dimacs/USA-road-d.COL.gr",
        "data/dimacs/USA-road-d.FLA.gr"
    };

    FILE *result_file=fopen("results/benchmark_results.csv","w");
    if(result_file==NULL){
        perror("results/benchmark_results.csv");
        return 1;
    }

    FILE *query_file=fopen("results/benchmark_queries.csv","w");
    if(query_file==NULL){
        perror("results/benchmark_queries.csv");
        fclose(result_file);
        return 1;
    }

    FILE *operation_file=fopen("results/benchmark_operations.csv","w");
    if(operation_file==NULL){
        perror("results/benchmark_operations.csv");
        fclose(result_file);
        fclose(query_file);
        return 1;
    }

    fprintf(result_file,
            "dataset,vertices,arcs,queries,binary_seconds,ternary_seconds,"
            "fibonacci_seconds,binomial_seconds,binary_avg_ms,"
            "ternary_avg_ms,fibonacci_avg_ms,binomial_avg_ms,"
            "matching_answers\n");
    fprintf(query_file,
            "dataset,query_id,source,target,binary_distance,"
            "ternary_distance,fibonacci_distance,binomial_distance\n");
    fprintf(operation_file,
            "dataset,heap,queries,scanned_arcs,settled_vertices,"
            "pq_inserts,pq_extracts,pq_decrease_keys\n");

    int ok=1;
    for(int i=0;i<3;i++){
        if(!benchmark_dataset(names[i],paths[i],i,
                              result_file,query_file,operation_file)){
            ok=0;
            break;
        }
    }

    fclose(result_file);
    fclose(query_file);
    fclose(operation_file);

    if(!ok){
        fprintf(stderr,"Benchmark failed or heap answers differed\n");
        return 1;
    }

    printf("Saved results/benchmark_results.csv, results/benchmark_queries.csv, "
           "and results/benchmark_operations.csv\n");
    return 0;
}
