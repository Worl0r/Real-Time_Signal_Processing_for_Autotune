#ifndef RING_BUFFER
#define RING_BUFFER

#include <stdlib.h>
#include <stdio.h>

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

class RingBuffer{
  public:
    RingBuffer(int size);
    void writeRingBuffer();
    int readRingBuffer();
    ~RingBuffer();

    int * dataRingBuffer;

  private:

    int sizeRingBuffer;
    int indexRead;
    int indexWrite;
};

class BufferDump{
    public:
        BufferDump();
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

#endif
