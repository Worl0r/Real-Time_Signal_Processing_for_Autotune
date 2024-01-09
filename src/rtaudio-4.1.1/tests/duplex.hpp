#ifndef RING_BUFFER
#define RING_BUFFER

#include <stdlib.h>
#include <stdio.h>

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

// Ring buffer for fundamental frequency
typedef struct {
  int sizeRingBuffer;
  int indexRead;
  int indexWrite;
  int* buffer;
} RingBuffer;

// Buffer to record some transformation through the test
typedef struct {
    MY_TYPE * bufferDump;
    int bufferFrameSize;
    int bufferDumpSize;
    int indexBufferDump;
    unsigned int bytes;
    std::string name;

} BufferOptions;

// Structure to store the autotune records
typedef struct{
  MY_TYPE * autoCorr;
  MY_TYPE * magnitudes;
  MY_TYPE * frequencies;
  MY_TYPE * oldFrequencies;
  MY_TYPE * phase;
  MY_TYPE * oldPhase;
  MY_TYPE * autotuneSignal;
} ProcessTools;

// Parametric structure for the autotune process
typedef struct {
  BufferOptions ** buffer;
  ProcessTools * processTools;
  RingBuffer * ringBufferFreq;
} Buffers;

// Functions
RingBuffer * createRingBuffer(int size);
int writeRingBuffer(RingBuffer * ringBuffer, int data);
void displayRingBuffer(RingBuffer * ringBuffer);
int getSizeRingBuffer(RingBuffer * ringBuffer);
int readRingBuffer(RingBuffer * ringBuffer);
void destroyRingBuffer(RingBuffer * ringBuffer);

#endif
