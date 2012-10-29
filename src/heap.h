#include <R.h>
#include <Rinternals.h>

typedef struct heap_element{ 
    int n;       /* dimension of E and F                */
    int *F;      /* set E                               */
    int *E;      /* set F                               */
    double b;    /* upper bound                         */
    int ibranch; /* branching index                     */
    int bound;   /* indicates "real" heap elements,     */
                 /* dummy elements have bound!=0        */
    int serial;  /* unique identifier for this element  */
    struct heap_element *up;   /* ptr to upper neigbour */
    struct heap_element *down; /* ptr to lower neigbour */
} heap_element_type;
               
typedef struct heap{
    struct heap_element *bottom; /* dummy element, lower entry point */
    struct heap_element *top;    /* dummy element, upper entry point */
    int count;                   /* # of elements                    */
    int maxcount;                /* max. # of elements               */
    int next_serial;             /* used to generate el->serial      */
} heap_type;

void int_copy(int *dest, int *src, int n);

heap_element_type *heap_element_init(int *F, int *E, int n,  double b, int ib);

void heap_element_copy(heap_element_type *dest, heap_element_type *orig);

void heap_element_destroy(heap_element_type *el);

void heap_element_print(heap_element_type *el);

heap_type *heap_init();

void heap_insert(heap_type *heap, heap_element_type *el);

void heap_traverse(heap_type *heap);

heap_element_type *heap_drop(heap_type *heap);

int heap_length(heap_type *heap);
