#include "printLock.h"

pthread_mutex_t printingLock = PTHREAD_MUTEX_INITIALIZER;

void printLock(){
    pthread_mutex_lock(&printingLock);
}

void printUnlock(){
    pthread_mutex_unlock(&printingLock);
}