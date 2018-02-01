// Created in conjunction with Paul Langton at Northeastern University

#include <stdint.h>
#include <sys/mman.h>
#include <assert.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>

#include "fast_malloc.h"
#define MIN(a,b) ((a) < (b) ? a : b)

typedef struct nu_free_cell {
    int64_t              size;
    struct nu_free_cell* next;
} nu_free_cell;


typedef struct bins {
    nu_free_cell* b32;
    nu_free_cell* b64;
    nu_free_cell* b128;
    nu_free_cell* b256;
    nu_free_cell* b512;
    nu_free_cell* b1024;
    nu_free_cell* b2048;
    nu_free_cell* b4096;
    nu_free_cell* b8192;
    nu_free_cell* b16384;
} mem_bins;

static const int64_t PAGE_SIZE = 4096;
static const int64_t CELL_SIZE  = (int64_t)sizeof(nu_free_cell);
static const int64_t MAX_BIN_SIZE = 16384;

__thread mem_bins* bins = NULL;
__thread int init = 0;

__thread pthread_mutex_t mutex;

void
bin_push(nu_free_cell** head, nu_free_cell* item) 
{
    item->next = *head;
    *head = item;
    assert((*head)->next != *head);
}

nu_free_cell*
bin_pop(nu_free_cell** head)
{
    nu_free_cell* mem = *head;
    *head = (*head)->next;
    assert(*head != mem);
    return mem;
}


int64_t
nu_free_list_length(nu_free_cell* n)
{
    int len = 0;

    for (nu_free_cell* pp = n; pp != 0; pp = pp->next) {
        len++;
    }

    return len;
}

void
nu_print_free_list(nu_free_cell* n)
{
    
    printf("Length: %d\n", nu_free_list_length(n));  
    for (; n != 0; n = n->next) {
        printf("%lx: (cell %ld %lx)\n", (int64_t) n, n->size, (int64_t) n->next); 

    }
}

void
print_bins()
{
    printf("= Bins: =\n");
    printf("32-Bytes:\n");
    nu_print_free_list(bins->b32);
    printf("64-Bytes:\n");
    nu_print_free_list(bins->b64);
    printf("128-Bytes:\n");
    nu_print_free_list(bins->b128);
    printf("256-Bytes:\n");
    nu_print_free_list(bins->b256);
    printf("512-Bytes:\n");
    nu_print_free_list(bins->b512);
    printf("1024-Bytes:\n");
    nu_print_free_list(bins->b1024);
    printf("2048-Bytes:\n");
    nu_print_free_list(bins->b2048);
    printf("4096-Bytes:\n");
    nu_print_free_list(bins->b4096);
    printf("8192-bytes:\n");
    nu_print_free_list(bins->b8192);
    printf("16384-bytes:\n");
    nu_print_free_list(bins->b16384);
    printf("= End of Bins =\n");
}


int64_t
get_bin_size(int64_t size) {
    return pow(2, ceil(log(size)/log(2)));
}

nu_free_cell**
get_bin(int64_t size)
{
    if (size <= 32) {
	    return &bins->b32;
    }
    else if (size <= 64) {
	    return &bins->b64;
    } 
    else if (size <= 128) {
	    return &bins->b128;
    }
    else if (size <= 256) {
	    return &bins->b256;
    }
    else if (size <= 512) {
        return &bins->b512;
    }
    else if (size <= 1024) {
	    return &bins->b1024;
    }
    else if(size <= 2048) {
	    return &bins->b2048;
    }
    else if(size <= 4096) {
	    return &bins->b4096;
    }
    else if (size <= 8192) {
	    return &bins->b8192;
    }
    else if (size <= 16384) {
        return &bins->b16384;
    }
    else {
	    return NULL;
    }
}



// Breaks down one block of the the given bin until the bin of the needed size has elements in
// it.
void
cascade(nu_free_cell** bin, int64_t size) 
{
    
    nu_free_cell* desired_bin = *get_bin(size);
    
    //print_bins();

    if (*bin == desired_bin) {
        return;
    }
    else {
        nu_free_cell* to_split = bin_pop(bin);
        int64_t split_size = to_split->size / 2;

        nu_free_cell* left_split = to_split;
        nu_free_cell* right_split = (nu_free_cell*) (((void*)to_split) + split_size);
        
        left_split->size = split_size;
        right_split->size = split_size;

        nu_free_cell** next_bin = get_bin(split_size);
        bin_push(next_bin, left_split);
        bin_push(next_bin, right_split);

        
        cascade(next_bin, size);
    }

}


// Gets the next biggest bin to split from. Fills the 4096 bin if nothing found
nu_free_cell**
get_next_biggest_bin(int64_t size) 
{
    int64_t bin_size = get_bin_size(size); 

    nu_free_cell** next_biggest_bin = NULL;
    for (int i = bin_size; i <= MAX_BIN_SIZE; i *= 2) {
        next_biggest_bin = get_bin(i);
        if (*next_biggest_bin != NULL) {
            break;
        }
    }

    if (*next_biggest_bin == NULL) {
        nu_free_cell** bin4096 = get_bin(MAX_BIN_SIZE);
        void* addr = mmap(0, MAX_BIN_SIZE, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANON, -1, 0);
        (*bin4096) = (nu_free_cell*) addr; 
        (*bin4096)->size = MAX_BIN_SIZE;
        (*bin4096)->next = NULL;
        return bin4096;
    }
    else {
        return next_biggest_bin;
    }
    
}


static
nu_free_cell*
get_cell(int64_t size)
{ 
    if (size > MAX_BIN_SIZE) {
        void* addr = mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        nu_free_cell* cell = (nu_free_cell*) addr; 
        cell->size = size;
        return cell;            
    }

    nu_free_cell** head = get_bin(size);
    if (*head == NULL) {
        cascade(get_next_biggest_bin(size), size);
    }

    return bin_pop(head);
}

static
nu_free_cell*
make_cell(int64_t size)
{
    void* addr = mmap(0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    nu_free_cell* cell = (nu_free_cell*) addr;
    cell->size = size;
    return cell;
}


void*
fast_malloc(size_t usize)
{
    pthread_mutex_init(&mutex, 0);
    pthread_mutex_lock(&mutex);

    int64_t size = (int64_t) usize;

    int64_t alloc_size = size + sizeof(int64_t);
    
    if (!init) {
        void* addr = mmap(0, sizeof(mem_bins), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        bins = (mem_bins*) addr;
        bins->b32 = NULL;
        bins->b64 = NULL;
        bins->b128 = NULL;
        bins->b256 = NULL;
        bins->b512 = NULL;
        bins->b1024 = NULL;
        bins->b2048 = NULL;
        bins->b4096 = NULL;
        bins->b8192 = NULL;
        bins->b16384 = NULL;


        for (int i = 0; i < 1; ++i) {
            nu_free_cell* page = make_cell(MAX_BIN_SIZE);
            bin_push(&bins->b16384, page);
        }
        init = 1;
    }

    printf("Mallocing %ld ...\n", usize);
    //print_bins();

   
    nu_free_cell* user_space = get_cell(alloc_size);
    
    *((int64_t*)user_space) = user_space->size;

    pthread_mutex_unlock(&mutex);
    return ((void*) user_space) + sizeof(int64_t);
}

void*
fast_realloc(void* prev, size_t bytes)
{    

    int64_t new_size = bytes;
    void* new_mem = fast_malloc(new_size);

    pthread_mutex_init(&mutex, 0);
    pthread_mutex_lock(&mutex);
    printf("Reallocating %p with %ld bytes...\n", prev, bytes);
    //print_bins();

    int64_t* prev_data = (int64_t*) prev;

    prev = prev - sizeof(int64_t);
    int64_t prev_size = *((int64_t*) prev);
    prev = prev + sizeof(int64_t);

    int64_t data_length = MIN(prev_size, new_size);

    memcpy(new_mem, prev, data_length);

    pthread_mutex_unlock(&mutex);

    fast_free(prev);

    return new_mem;
 
}

void
fast_free(void* addr) 
{    
    printf("Freeing %p...\n", addr);
    //print_bins();
    nu_free_cell* cell = (nu_free_cell*)(addr - sizeof(int64_t));
    int64_t size = *((int64_t*) cell);

    if (size > MAX_BIN_SIZE) {
        munmap((void*) cell, size);
    }
    else {
        pthread_mutex_init(&mutex, 0);
        pthread_mutex_lock(&mutex);
  
        cell->size = size;
        nu_free_cell** head = get_bin(size);
        bin_push(head, cell);


        pthread_mutex_unlock(&mutex);
    }

}

