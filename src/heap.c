/* #include "stdlib.h" */
/* #include "math.h" */

#include "heap.h"
/*
 * #define N 10
 * #define M 20
 */


void int_copy(int *dest, int *src, int n)
{
    int i;
    for(i=0;i<=n-1;i++)
	dest[i]=src[i];
}

/* 
 * element functions: 
 * init
 * copy 
 * destroy
 * print
 */
heap_element_type *heap_element_init(int *F, int *E, int n,  double b, int ib)
{
    heap_element_type *heap_element;

    /* allocate memory, initialize values */
    heap_element = Calloc((size_t)1,heap_element_type);
    if(heap_element == NULL )
	{
	    warning("Calloc failed for heap element");
	    return(NULL);
	};
    heap_element->n      = n;
    heap_element->F = Calloc((size_t)n,int);
    if(heap_element->F == NULL)
	{
	    warning("Calloc failed for F");
	    return(NULL);
	};
    heap_element->E = Calloc((size_t)n,int);
    if(heap_element->E == NULL )
	{
	    warning("Calloc failed for E");
	    return(NULL);
	};
    if(F!=NULL) int_copy(heap_element->F,F,n);
    if(E!=NULL) int_copy(heap_element->E,E,n);
    heap_element->b      = b;
    heap_element->bound  = 0;
    heap_element->ibranch=ib;
    heap_element->serial = -1; /* -1 indicates a initialised only element */
    heap_element->up     = NULL;
    heap_element->down   = NULL;
    return(heap_element);
}

void heap_element_copy(heap_element_type *dest, heap_element_type *orig)
{
    int_copy(dest->F, orig->F, orig->n);
    int_copy(dest->E, orig->E, orig->n);
    dest->n       = orig->n;
    dest->b       = orig->b;
    dest->bound   = orig->bound;
    dest->ibranch = orig->ibranch;
    dest->up      = orig->up;
    dest->down    = orig->down;    
    dest->serial  = -2; /* -2 indicates a copied element */
}

void heap_element_destroy(heap_element_type *el)
{
    Free(el->F);
    Free(el->E);
    Free(el);
}

void heap_element_print(heap_element_type *el)
{
#ifdef HEAP_DEBUG
    int i,up=-1,down=-1;
    printf("element %i, n=%i, bound=%f\nF: ",el->serial, el->n, el->b);
    for(i=0;i<=(el->n)-1;i++)
	printf("%2i ",el->F[i]);
    printf("\nE: ");
    for(i=0;i<=(el->n)-1;i++)
	printf("%2i ",el->E[i]);
    printf("\n");
    if(el->up)   up=el->up->serial;
    if(el->down) down=el->down->serial;
    printf("pointers: up->%i, down->%i\n",up,down);
#endif
}

/* 
 * heap funcions:
 * init
 * insert
 * traverse
 * drop
 * length
 */
heap_type *heap_init()
{
    heap_type *heap;
	int *i,*j,*k,*l;
    /* dummy arrays */
	i=Calloc(1,int);
	j=Calloc(1,int);
	k=Calloc(1,int);
	l=Calloc(1,int);
    i[1]=0;
    j[1]=0;
    k[1]=0;
    l[1]=0;

    /* 
     * allocate memory, initialize values 
     */
    heap = Calloc((size_t)1,heap_type);
    if(heap == NULL)
	{
	    warning("Calloc failed for heap"); 
	    return(NULL);
	}

    /* 
     *  lower and upper entry to the heap: 
     */
    heap->bottom = heap_element_init(i,j,1,0,0);
    heap->top    = heap_element_init(k,l,1,0,0);

    /* 
     * initialize links between elements 
     */
    heap->bottom->up = heap->top;
    heap->top->down  = heap->bottom;  

    /* 
     * indicate top and bottom  by bound=+/-1 
     */
    heap->top->bound     =  1;
    heap->bottom->bound  = -1;
    heap->count          =  0;
    heap->maxcount       =  0;
    heap->top->serial    =  0; /* 0 for boundary elements */
    heap->bottom->serial =  0; /*                         */
    heap->next_serial    =  1;
    return(heap);
}

void heap_insert(heap_type *heap, heap_element_type *el)
{
    heap_element_type *ptr;
    
    ptr=heap->bottom->up;

    /* search from bottom to top where el fits in the heap */
    if (ptr->bound == +1)
	/* insert first element between top and bottom */	
	{
#ifdef HEAP_DEBUG
	    /* printf("inserting first element with bound %f\n",el->b); */
#endif
	    /* update new links from/to el */
	    el->up        = ptr;
	    el->down      = ptr->down;	     
	    ptr->down->up = el;
	    ptr->down     = el;
	}
    else
	{
#ifdef HEAP_DEBUG
	    /* printf("inserting element with bound %f\n",el->b); */
#endif
	    if (el->b <= ptr->b)		
		/* insert below lowermost element */
		{
		    /* update new links from/to el */
		    el->up        = ptr;
		    el->down      = ptr->down;	     
		    ptr->down->up = el;
		    ptr->down     = el;
		}
	    else
		{
		    /* go up until top or bound el->b */
		    while(ptr->bound != +1 && ptr->b <= el->b)
			ptr=ptr->up;
		    /* insert below found element */
		    /* update new links from/to el */

		    el->up           = ptr;
		    el->down         = ptr->down;	    
		    ptr->down->up    = el;
		    ptr->down        = el;
		}
	}
    el->serial        = heap->next_serial;
    heap->next_serial = heap->next_serial+1;;
    heap->count       = heap->count+1;
    if( heap->count > heap->maxcount) heap->maxcount = heap->count;
}

void heap_traverse(heap_type *heap)
{
    heap_element_type *ptr;
    int i=1;

    ptr=heap->bottom->up;
#ifdef HEAP_DEBUG
	    printf("-------- heap contents: --------\n",i);
#endif

    while(ptr->bound != +1 && i < 20)
	{
#ifdef HEAP_DEBUG
	    printf("%i th heap element:\n",i);
#endif
	    heap_element_print(ptr);
	    ptr=ptr->up;
	    i++;
	}
#ifdef HEAP_DEBUG
	    printf("--------------------------------\n",i);
#endif
}

heap_element_type *heap_drop(heap_type *heap)
{
    heap_element_type *ptr,*destroy;
    int lerr;

    if(! (ptr=heap_element_init(NULL,NULL,heap->top->down->n,0,0))){
	/* error("memory exhausted"); exit(1); */
        lerr=1;
    }

    if(heap->top->down->bound!=-1)
	{
	    heap_element_copy(ptr,heap->top->down);
	    
	    destroy                   = heap->top->down;
	    heap->top->down->down->up = heap->top;
	    heap->top->down           = heap->top->down->down;

	    heap_element_destroy(destroy);
	    heap->count = heap->count-1;
	    return(ptr);
	}
    else
	return(NULL);
}

int heap_length(heap_type *heap)
{
    return(heap->count);
}

/* cc -I/usr/lib/R/include/ heap.c -lm -g -DHEAP_DEBUG=1 \
 *    -L/usr/lib/R/bin -lR 
 * LD_LIBRARY_PATH=/usr/lib/R/bin ./a.out  
 */

/*
int M=200;
int N=40;

 int main(int argc, void *argv) 
{
    
        heap_type *heap; 
        heap_element_type *el; 
        int i; 
        int f[M]; 
        int e[M]; 
   
        heap = heap_init();     
   
        for(i=1;i<=M;i++) 
     {
	 
	            f[i]=-1; 
	            e[i]=i; 
     }
    
        for(i=1;i<=N;i++) 
     {
	 
	            el = heap_element_init(f,e,M,rand()*(double)sin((double)i),i); 
	            heap_insert(heap, el); 
     }
    
   
        heap_traverse(heap); 
   
        for(i=N;i>=-111;i--) 
     {
	 
	 #ifdef HEAP_DEBUG 
	            printf("removing number %i:\n",i); 
	 #endif 
	            el=heap_drop(heap); 
	            if(el!=NULL) 
	                  heap_element_print(el); 
	            else 
	   #ifdef HEAP_DEBUG 
	                  printf("element is NIL!\n"); 
	 #endif 
	            heap_traverse(heap); 
     }
    
   
        exit(0);
}
 
*/
     
