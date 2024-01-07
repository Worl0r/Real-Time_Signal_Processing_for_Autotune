#ifndef RING_BUFFER
#define RING_BUFFER

#include <stdlib.h>
#include <stdio.h>

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

// Some useful functions
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// Functions

int fft(double *x, double *y, const int m);
char *getmem(int leng, unsigned size);
double *dgetmem(int leng);
static int checkm(const int m);

// Struct
// class RingBuffer{
//   public:
//     RingBuffer(int size);
//     int writeRingBuffer(int data);
//     void displayRingBuffer();
//     int getSizeRingBuffer();
//     int readRingBuffer();
//     ~RingBuffer();

//   private:
//     int sizeRingBuffer;
//     int indexRead;
//     int indexWrite;
//     int * buffer;
// };


typedef struct {
  int sizeRingBuffer;
  int indexRead;
  int indexWrite;
  int* buffer;
} RingBuffer;

typedef struct {
    MY_TYPE * bufferDump;
    int bufferFrameSize;
    int bufferDumpSize;
    int indexBufferDump;
    unsigned int bytes;
    std::string name;

} BufferOptions;

typedef struct{
  MY_TYPE * auto_corr;
  MY_TYPE * im_input;
  MY_TYPE * magnitudes;
  MY_TYPE * frequencies;
  MY_TYPE * phase;
  MY_TYPE * oldPhase;
  MY_TYPE * autotuneSignal;
} ProcessTools;

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
