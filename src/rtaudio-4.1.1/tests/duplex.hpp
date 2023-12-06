#ifndef RING_BUFFER
#define RING_BUFFER

#include <stdlib.h>
#include <stdio.h>

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

class RingBuffer{
  public:
    RingBuffer(int size);
    int writeRingBuffer(int data);
    void displayRingBuffer();
    MY_TYPE readRingBuffer();
    ~RingBuffer();

  private:
    int sizeRingBuffer;
    int indexRead;
    int indexWrite;
    MY_TYPE * buffer;
};

typedef struct {
    MY_TYPE * bufferDump;
    int bufferFrameSize;
    int bufferDumpSize;
    int sampleFrenquency;
    int indexBufferDump;
    unsigned int bytes;
    std::string name;

} BufferOptions;

typedef struct {
  BufferOptions *buffer;
} Buffers;

#endif
