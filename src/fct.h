#include <stdlib.h>
#include <stdio.h>

// Stucture
typedef int * listRingBuffer;

typedef struct {
    listRingBuffer list;
    int indexWrite;
    int indexRead;
    int sizeRingBuffer;
}ringBuffer ;

// Les fonctions
void printRingBuffer( ringBuffer ringBufferData);
void addRingBuffer(int data, ringBuffer * ptr_ringBufferData);
int getRingBuffer(ringBuffer * ptr_ringBufferData);
ringBuffer * createRingBuffer(int n);
