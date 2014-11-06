/*
 * heap.c
 *
 * Implements a static heap, i.e. the allocation of the heap is fixed
 * and set at creation.  See <heap.h> for utility functions.
 *
 * Implementation derived from Chapter 7: Heapsort of
 * "Introdution to Algorithms" by Cormen, Leiserson and Rivest,
 * The MIT Press.
 */

#include <stdio.h>
#include "mex.h"

#include "cityheap.h"

#define PARENT(i)	(((i)-1) >> 1)
#define LEFT(i)		(((i) << 1) + 1)
#define RIGHT(i)	(((i)+1) << 1)

void build_heap(Heap *heap);
void heapify(Heap *heap, int i);


/*
 * Creates an empty array for an element type for which greater_than()
 * is defined over.
 */
Heap *make_empty_heap(initial_heap_size)
int	initial_heap_size;
{
  Heap	*heap;

  if (!(heap = (Heap *) mxCalloc( 1, sizeof(Heap) ))) {
    mexErrMsgTxt("make_empty_heap(): unable to allocate memory for heap structure.\n");
  }
     
  heap->heap_max_size = initial_heap_size;
  heap->heap_size = 0;
     
  if (!(heap->array = (active_elem **) mxCalloc(initial_heap_size, sizeof(active_elem *)))) {
    mexErrMsgTxt("make_empty_heap(): unable to allocate memory for heap array.\n");
  }

  return( heap );
}

/*
 * Creates a heap from a given array.
 */
Heap *make_heap(array, size)
     active_elem    **array;
     int	size;
{
     Heap	*heap;

     if (!(heap = (Heap *) mxCalloc(1, sizeof(Heap) ))) {
	  mexErrMsgTxt("make_heap(): unable to allocate memory for heap structure.\n");
     }

     heap->heap_max_size = size;
     heap->heap_size = size;
     heap->array = array;			/* use the same array */

     /* build the heap */
     build_heap(heap);

     return( heap );
}


/*
 * Free dynamically allocated structures in heap.
 */
void free_heap(heap)
     Heap	*heap;
{
     if (heap && heap->array) mxFree(heap->array);
     mxFree(heap);
}


/*
 * Assumes left and right subtree of root are heaps only that the root node
 * may not be the largest element in the heap and therefore needs to be
 * "sunked" to the right position in the heap.  Runs in O(lg n).
 */
void heapify(heap, i)
     Heap		*heap;
     register int	i;
{
     register int	l, r, largest;
     register int	heap_size = heap->heap_size;
     register active_elem	**heap_array = heap->array, *tmp;

     for (;;) {

	  l = LEFT(i); r = RIGHT(i);
	  largest = ((l < heap_size) && (heap_array[l]->dist < heap_array[i]->dist)) ? l : i;
	  if ((r < heap_size) && (heap_array[r]->dist < heap_array[largest]->dist)) largest = r;
	  if (largest == i) break;
	  else {
	       tmp = heap_array[i]; heap_array[i] = heap_array[largest]; heap_array[largest] = tmp;
	       i = largest;
	  }
     }
}

/*
 * Takes an ordinary array and makes it into a heap.  Runs in O(n).
 */
void build_heap(heap)
     Heap	*heap;
{
     register int	i;
     
     for (i=(heap->heap_size >> 1) - 1; i>=0; i--)
	  heapify(heap, i);
}

/*
 * Returns a pointer to the maximal element in the heap
 * and removes it from the heap as well.  Returns NULL
 * if there are no more elements in the heap.
 * Runs in O(lg n).
 */

active_elem *heap_extract_max(heap)
     Heap	*heap;
{
     register active_elem	**heap_array = heap->array, *max;

     if (heap->heap_size <= 0) return( NULL );
     max = heap_array[0];
     heap_array[0] = heap_array[--heap->heap_size];
     heapify(heap, 0);

     return( max );
}

/*
 * Inserts the new element into the heap preserving
 * all properties of a heap.  Runs in O(lg n).
 */
void heap_insert(heap, elt)
     Heap	*heap;
     active_elem	*elt;
{
     register int	i, p;
     register active_elem	**heap_array = heap->array;

     if (heap->heap_size >= heap->heap_max_size) {
	  mexErrMsgTxt("heap_insert():  heap overflow.\n");
     }
     
     i = heap->heap_size++;
     while ((i > 0) && (elt->dist < heap_array[p=PARENT(i)]->dist)) {
	  heap_array[i] = heap_array[p];
	  i = p;
     }
     
     heap_array[i] = elt;
}

