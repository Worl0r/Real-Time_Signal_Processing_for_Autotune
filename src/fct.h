#include <stdlib.h>
#include <stdio.h>

// Types
typedef int * listRingBuffer;

typedef struct {
    listRingBuffer list;
    int indexWrite;
    int indexRead;
    int sizeRingBuffer;
}ringBuffer ;

// Functions
void printRingBuffer( ringBuffer ringBufferData);
void writeRingBuffer(int data, ringBuffer * ptr_ringBufferData);
int readRingBuffer(ringBuffer * ptr_ringBufferData);
ringBuffer * createRingBuffer(int n);
