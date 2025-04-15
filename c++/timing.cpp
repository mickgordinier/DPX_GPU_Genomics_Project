#include "timing.h"

struct timeval tv;
uint64_t start_time;

uint64_t start_timer(){
    gettimeofday(&tv,NULL);
    start_time = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    return start_time;
}

uint64_t get_time(){
    gettimeofday(&tv,NULL);
    uint64_t current_time = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    return current_time;
}

uint64_t get_elapsed_time(){
    gettimeofday(&tv,NULL);
    uint64_t current_time = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    uint64_t elapsed_time = current_time - start_time;
    return elapsed_time;
}