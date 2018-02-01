#ifndef FASTMALLOC_H
#define FASTMALLOC_H

void* fast_malloc(size_t usize);
void* fast_realloc(void* prev, size_t bytes);
void fast_free(void* addr);

#endif
