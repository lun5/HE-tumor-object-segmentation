/*
 * heap.h
 */

/* 
 * Heap data structure
 */

typedef struct {
     int coord[3];
     float dist;
     int prevdir;
   } active_elem;

typedef struct {
     active_elem    **array;				/* heap array */
     int	heap_max_size;				/* size of heap array */
     int	heap_size;				/* number of entries in array used */
} Heap;


/* 
 * Heap utility functions
 */
extern Heap 	*make_empty_heap(int initial_heap_size);
extern Heap 	*make_heap(active_elem **array, int size);
#define 	heap_empty_p(heap)	((heap)->heap_size <= 0)
#define 	heap_full_p(heap)	((heap)->heap_size >= (heap)->heap_max_size)
extern void 	free_heap(Heap *heap);

extern void 	heap_insert(Heap *heap, active_elem *elt);
extern active_elem 	*heap_extract_max(Heap *heap);
#define 	heap_peek_max(heap)	((heap)->array[0])
