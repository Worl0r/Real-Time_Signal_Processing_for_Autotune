#include <stdlib.h>
#include <stdio.h>
#include "fct.h"

ringBuffer * createRingBuffer(int n)
{
    ringBuffer * ptr_ringBuffer =  malloc(sizeof(ringBuffer));
    ptr_ringBuffer->list = malloc(sizeof(int)*n);
    ptr_ringBuffer->indexWrite = 0;
    ptr_ringBuffer->indexRead = -1;
    ptr_ringBuffer->sizeRingBuffer = n;

    return ptr_ringBuffer;
}

void writeRingBuffer(int data, ringBuffer * ptr_ringBufferData){
    if (ptr_ringBufferData->indexWrite == ptr_ringBufferData->indexRead ){
        ptr_ringBufferData->indexRead++;
        ptr_ringBufferData->indexRead = ptr_ringBufferData->indexRead % ptr_ringBufferData->sizeRingBuffer;
    }

    if (ptr_ringBufferData->indexRead == -1){
        ptr_ringBufferData->indexRead = 0;
    }

    ptr_ringBufferData->list[ptr_ringBufferData->indexWrite] = data;
    ptr_ringBufferData->indexWrite++;
    ptr_ringBufferData->indexWrite = ptr_ringBufferData->indexWrite % ptr_ringBufferData->sizeRingBuffer;
}

void printRingBuffer(ringBuffer ringBufferData){
    int n = ringBufferData.sizeRingBuffer;

    printf("Ring Buffer is :\n");
    for(int i=0; i <n; i++){
        printf("%i\n",ringBufferData.list[i]);
    }

    printf("indexRead is: %i\n",ringBufferData.indexRead);
    printf("indexWrite is: %i\n",ringBufferData.indexWrite);

}

int readRingBuffer(ringBuffer * ptr_ringBufferData){

    if (ptr_ringBufferData->indexRead == -1){
        return -1;
    }

    int value = ptr_ringBufferData->list[ptr_ringBufferData->indexRead];
    ptr_ringBufferData->indexRead++;
    ptr_ringBufferData->indexRead = ptr_ringBufferData->indexRead % ptr_ringBufferData->sizeRingBuffer;
    return value;
}
