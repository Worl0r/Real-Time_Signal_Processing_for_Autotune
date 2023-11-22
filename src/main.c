#include <stdlib.h>
#include <stdio.h>
#include "fct.c"


int main(void)
{
    int sizeRingBuffer = 5;
    ringBuffer * ptr_ringBufferData =  createRingBuffer(sizeRingBuffer);

    // Test0
    printRingBuffer(*ptr_ringBufferData);

    // Test1
    addRingBuffer(1,  ptr_ringBufferData);
    printRingBuffer(*ptr_ringBufferData);

    // Test2
    addRingBuffer(2,  ptr_ringBufferData);
    addRingBuffer(3,  ptr_ringBufferData);
    addRingBuffer(4,  ptr_ringBufferData);
    addRingBuffer(5,  ptr_ringBufferData);
    addRingBuffer(6,  ptr_ringBufferData);
    printRingBuffer(*ptr_ringBufferData);

    // Test3
    addRingBuffer(7,  ptr_ringBufferData);
    printRingBuffer(*ptr_ringBufferData);

    // Test4
    int getItemRingBuffer = getRingBuffer( ptr_ringBufferData);
    printRingBuffer(*ptr_ringBufferData);
    printf("L element est : %i\n", getItemRingBuffer);

    // Test5
    ringBuffer * ptr_ringBufferData_test5 =  createRingBuffer(sizeRingBuffer);
    int getItemRingBuffer_test5 = getRingBuffer( ptr_ringBufferData_test5);
    printRingBuffer(*ptr_ringBufferData_test5);
    printf("getItemRingBuffer : %i\n",getItemRingBuffer_test5);

    // Test6
    addRingBuffer(1,  ptr_ringBufferData_test5);
    printRingBuffer(*ptr_ringBufferData_test5);
    getItemRingBuffer_test5 = getRingBuffer( ptr_ringBufferData_test5);
    printf("getItemRingBuffer : %i\n",getItemRingBuffer_test5);
    printRingBuffer(*ptr_ringBufferData_test5);

    return 1;
}