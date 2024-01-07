/******************************************/
/*
  duplex.cpp
  by Gary P. Scavone, 2006-2007.

  This program opens a duplex stream and passes
  input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <math.h>
#include <ctime>
#define _USE_MATH_DEFINES
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>

#include "duplex.hpp"

/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8
*/

// typedef signed short MY_TYPE;
// #define FORMAT RTAUDIO_SINT16

using namespace std;

/*
typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32
*/

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device to use (default = 0),\n";
  std::cout << "    oDevice = optional output device to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
  exit( 0 );
}

/////////////////////////////////////////////// Options ///////////////////////////////////////////////

const string PATH_RECORD = "../../../output/";
const int bufferSize = 100000;
const int ringBufferSize = 20;
double *_sintbl = 0;
int maxfftsize = 512;
int samplingFrequency = 48000;
int nbrHarmonics = 5;
int nbrDemiTons = 2;

/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////// Utils ///////////////////////////////////////////////


RingBuffer* createRingBuffer(int size) {
  RingBuffer* rb = (RingBuffer*) malloc(sizeof(RingBuffer));
  rb->sizeRingBuffer = size;
  rb->indexRead = -1;
  rb->indexWrite = 0;
  rb->buffer = (int*) calloc(size, sizeof(int));
  return rb;
}

void destroyRingBuffer(RingBuffer* rb) {
  free(rb->buffer);
  free(rb);
}

int writeRingBuffer(RingBuffer* rb, int data) {
  if (rb->indexRead == -1){
    rb->indexRead = 0;
  }

  rb->buffer[rb->indexWrite] = data;
  rb->indexWrite++;
  rb->indexWrite = rb->indexWrite % rb->sizeRingBuffer;

  return 0;
}

int readRingBuffer(RingBuffer* rb) {
  if (rb->indexRead == -1){
    printf("[ERROR] Buffer is empty in readRingBuffer\n");
    return -1;
  }

  if (rb->indexRead == rb->indexWrite){
    printf("[ERROR] indexRead in ring buffer overflow the buffer\n");
    return -1;
  }

  int data = rb->buffer[rb->indexRead];
  rb->indexRead++;
  rb->indexRead = rb->indexRead % rb->sizeRingBuffer;
  return data;
}

int meanRingBuffer(RingBuffer* rb) {
  if (rb->indexRead == -1){
    printf("[ERROR] Buffer is empty in readRingBuffer\n");
    return -1;
  }

  int sum = 0;
  int coef = 0;

  for (int i = 0; i < rb->sizeRingBuffer; i++){
    if(rb->buffer[i] != 0){
      coef++;
      sum = sum + rb->buffer[i];
    }
  }
  rb->indexRead++;
  rb->indexRead = rb->indexRead % rb->sizeRingBuffer;

  return sum / coef;
}

int getSizeRingBuffer(RingBuffer* rb) {
  return rb->sizeRingBuffer;
}

void displayRingBuffer(RingBuffer* rb) {
  if (rb->indexRead == -1){
    printf("[ERROR] Buffer is empty in displayRingBuffer\n");
    return;
  }

  time_t result = time(NULL);
  printf("[LOG][TEST] You display your ring buffer at %s\n", asctime(localtime(&result)));
  printf("[ ");

  for(int i=0; i < rb->sizeRingBuffer; i++){
    if (i == rb->sizeRingBuffer-1){
      printf("%d", rb->buffer[i]);
    }
    else{
      printf("%d | ", rb->buffer[i]);
    }
  }

  printf(" ]\n");
  printf("Read index is : %d\n", rb->indexRead);
  printf("Write index is : %d\n", rb->indexWrite);
}

void displayTab(string testName, MY_TYPE * buffer, int size){
  cout << "[LOG] [TEST] " << testName << endl << "[";

  for (int i = 0; i<size; i++){
    cout << (double) buffer[i] << " | ";
  }
  cout << " ]" << endl;
}

int write_buff_dump(double* buff, const int n_buff, double* buff_dump, const int n_buff_dump, int* ind_dump) {

  int i = 0;
  for (i = 0; i < n_buff; i++)
  {
    if (*ind_dump < n_buff_dump)
    {
      buff_dump[*ind_dump] = buff[i];
      (*ind_dump)++;
    } else
    {
      break;
    }
  }

  return i;
}

void deallocateBuffer(BufferOptions * bufferOptions){
  free(bufferOptions->bufferDump);
  free(bufferOptions);
};

BufferOptions * allocateBuffer(int bufferDumpSize,
  int bufferFrameSize,
  unsigned int bufferBytes,
  string name,
  int indexBufferDump = 0){

  BufferOptions * bufferOptions = (BufferOptions *) malloc(sizeof(BufferOptions));
  if (bufferOptions == NULL) {
    fprintf(stderr, "Memory allocation error for BufferOptions\n");
    exit(1);
  }

  bufferOptions->bufferDumpSize = bufferDumpSize;

  bufferOptions->bufferDump = (MY_TYPE*) calloc(bufferOptions->bufferDumpSize, sizeof(MY_TYPE));
  if (bufferOptions->bufferDump == NULL) {
    fprintf(stderr, "Memory allocation error for bufferDump\n");
    exit(1);
  }

  bufferOptions->bufferFrameSize = bufferFrameSize;
  bufferOptions->indexBufferDump = indexBufferDump;
  bufferOptions->bytes = bufferBytes;
  bufferOptions->name =  name;

  return bufferOptions;
};

void writeBuffer(BufferOptions *bufferIn, string PATH_RECORD){

  cout << "[INFO] Saving of Signals: " << (bufferIn->name).c_str() << endl;
  FILE* f;
  f = fopen(strcat((char *) PATH_RECORD.c_str(), (char *) bufferIn->name.c_str()), "wb");

  assert(f);

  fwrite(bufferIn->bufferDump, sizeof(*(bufferIn->bufferDump)), bufferIn->bufferDumpSize, f);
  fclose(f);
};

// The goal is to extract the second maximum for the autocor function.
int extractFundamentalFrequency(MY_TYPE * autocor, int sizeAutocor){
  assert(autocor);
  assert(sizeAutocor >= 5);

  int * maxList = (int *) malloc(sizeof(int) * (sizeAutocor));
  int count = 0;

  // TODO: adapt this hard coded variable
  for (int i = 4; i < sizeAutocor - 1; i++){
    if ((autocor[i-1] < autocor[i]) && (autocor[i] > autocor[i+1])){
      maxList[count] = i,
      count++;
    }
  }

  int max = -1;

  for(int i = 0; i< sizeAutocor; i++){
    if (autocor[maxList[i]] > max){
      max = maxList[i];
    }
  }

  free(maxList);

  return max;
}

void demi_auto_corr(double * input, int size, MY_TYPE * auto_corr) {
  int n;
  int k;

  for (n = 0; n < size; n++) {
    auto_corr[n] = 0;
    for (k = n; k < size; k++) {
      auto_corr[n] += input[k] * input[k - n];
    }
  }
}

void hzTodemiTon(MY_TYPE * frequencies, int size){
  for (int i = 0; i < size; i++){
    frequencies[i] = 12 * log2(frequencies[i]);
  }
}

void demiTonToHz(MY_TYPE * frequencies, int size){
  for (int i = 0; i < size; i++){
    frequencies[i] = pow(2, frequencies[i]/12);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////// process //////////////////////////////////////////////
void changeFundamentalFrequency(MY_TYPE * fundamentalFrequency, int nbrDemiTons, int size){
  hzTodemiTon(fundamentalFrequency, size);

  //*fundamentalFrequency = round(*fundamentalFrequency / nbrDemiTons) * nbrDemiTons;

  demiTonToHz(fundamentalFrequency, size);
}

int saveFundamentalFrequency(RingBuffer * ringBuffer, ProcessTools * processTools, double * input, int size){

  // Clear autocorrelation table
  memset(processTools->auto_corr, 0, size * sizeof(MY_TYPE));

  // Compute autocorrelation
  demi_auto_corr(input, size, processTools->auto_corr);

  int index = extractFundamentalFrequency(processTools->auto_corr, size);

  if (index == -1){
    fprintf(stderr, "Fundamental frequency not found\n");
    exit(1);
  }

  // Save fundamental frequency
  writeRingBuffer(ringBuffer, index);

  return 1;
}

void process(double * input, RingBuffer * ringBufferFreq, ProcessTools * processTools, int size) {

  memset(processTools->im_input, 0, size * sizeof(MY_TYPE));

  // Compute fundamental frequency
  int val = saveFundamentalFrequency(ringBufferFreq, processTools, input, size);

  if (val == -1){
    fprintf(stderr, "Fundamental frequency not found\n");
    exit(1);
  }

  // FFT transformation
  int success = fft(input, processTools->im_input, size);

  if(success == -1) {
    fprintf(stderr, "FFT failed\n");
    exit(1);
  }

  // Compute mean of fundamental frequency
  int meanFreq = meanRingBuffer(ringBufferFreq);

  if(meanFreq == 0){
    return;
  }

  // Convert meanFreq to Hz
  MY_TYPE meanFreqHz = (MY_TYPE) (1.0 / meanFreq) * samplingFrequency;

  // Display
  displayRingBuffer(ringBufferFreq);
  cout << "frequencies_value[0] in Hz before T" << meanFreqHz << endl;
  cout << "frequencies_value[0] in idx before T" << meanFreq << endl;

  // Transform fundamental frequency
  changeFundamentalFrequency(&meanFreqHz, nbrDemiTons, nbrHarmonics);

  // Setup fundamental frequency and magnitude
  processTools->frequencies[0] = meanFreqHz;
  meanFreq = round((1.0 / meanFreqHz) * samplingFrequency);
  processTools->magnitudes[0] = (MY_TYPE) sqrt((input[meanFreq] * input[meanFreq]) + (processTools->im_input[meanFreq] * processTools->im_input[meanFreq]));
  //processTools->phase[0] = atan2(processTools->im_input[meanFreq], input[meanFreq]);
  processTools->phase[0] = processTools->oldPhase[0] + 2.0 *  M_PI * processTools->frequencies[0];

  cout << "frequencies_value[0] in Hz after T" << meanFreqHz << endl;
  cout << "frequencies_value[0] in idx after T" << meanFreq << endl;

  // Idem for other harmonics
  for (int j = 2; j < nbrHarmonics; j++){
    processTools->frequencies[j-1] = j * meanFreqHz;
    processTools->magnitudes[j-1] = sqrt((input[j * meanFreq] * input[j * meanFreq]) + (processTools->im_input[j * meanFreq] * processTools->im_input[j * meanFreq]));
    processTools->phase[j-1] = processTools->oldPhase[j-1] + 2.0 *  M_PI * processTools->frequencies[j-1];
    //processTools->phase[j-1] = atan2(processTools->im_input[j * meanFreq], input[j * meanFreq]);
  }

  displayTab("Magnitudes", processTools->magnitudes, nbrHarmonics);
  displayTab("Frequencies", processTools->frequencies, nbrHarmonics);
  displayTab("oldPhase", processTools->oldPhase, nbrHarmonics);
  displayTab("Phase", processTools->phase, nbrHarmonics);

  // Update old phase
  memcpy(processTools->oldPhase, processTools->phase, nbrHarmonics * sizeof(MY_TYPE));

  for(int i = 0; i< size; i++){
    MY_TYPE sum = 0;

    for(int j = 0; j < nbrHarmonics; j++){
      sum = sum + 2 * processTools->magnitudes[j] * cos(2.0 * M_PI * processTools->frequencies[j] * i + processTools->phase[j]);
    }

    if(sum > 1){
      processTools->autotuneSignal[i] = 0;
    }
    else if(sum < -1){
      processTools->autotuneSignal[i] = 0;
    }
    else{
      processTools->autotuneSignal[i] = sum;
    }
  }

  displayTab("SignalTab", (processTools->autotuneSignal), size);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  Buffers * listBuffers = (Buffers *) data;

  BufferOptions * bufferStructIn = listBuffers->buffer[0];
  BufferOptions * bufferStructOut = listBuffers->buffer[1];
  BufferOptions * autocor = listBuffers->buffer[2];
  BufferOptions * fundamentalFrequency = listBuffers->buffer[3];

  RingBuffer * ringBufferFreq = listBuffers->ringBufferFreq;

  ProcessTools * processTools = listBuffers->processTools;

  memcpy(outputBuffer, inputBuffer, (size_t) (bufferStructIn->bytes));

  // Record input buffer
  int index = write_buff_dump(
    (double *) inputBuffer,
    bufferStructIn->bufferFrameSize,
    (MY_TYPE *) bufferStructIn->bufferDump,
    bufferStructIn->bufferDumpSize,
    &bufferStructIn->indexBufferDump
  );

  /////////////////// Process ////////////////////////

  process((double *) inputBuffer, ringBufferFreq, processTools, (int) bufferStructIn->bufferFrameSize);

  ////////////////////////////////////////////////////

  // Record output buffer
  int index2 = write_buff_dump(
    (double *) (processTools->autotuneSignal),
    bufferStructIn->bufferFrameSize,
    (MY_TYPE *) bufferStructOut->bufferDump,
    bufferStructOut->bufferDumpSize,
    &bufferStructOut->indexBufferDump
  );

  return 0;
}

int main( int argc, char *argv[] )
{

  unsigned int channels, fs, bufferBytes, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

  // Minimal command-line checking
  if (argc < 3 || argc > 7 ) usage();

  RtAudio adac;
  if ( adac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  if ( argc > 3 )
    iDevice = (unsigned int) atoi(argv[3]);
  if ( argc > 4 )
    oDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    iOffset = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    oOffset = (unsigned int) atoi(argv[6]);

  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 512;
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = iDevice;
  iParams.nChannels = channels;
  iParams.firstChannel = iOffset;
  oParams.deviceId = oDevice;
  oParams.nChannels = channels;
  oParams.firstChannel = oOffset;

  if ( iDevice == 0 )
    iParams.deviceId = adac.getDefaultInputDevice();
  if ( oDevice == 0 )
    oParams.deviceId = adac.getDefaultOutputDevice();

  RtAudio::StreamOptions options;
  //options.flags |= RTAUDIO_NONINTERLEAVED;

  bufferBytes = bufferFrames * channels * sizeof( MY_TYPE );

  /////////////////////////////////////////////// Buffers ////////////////////////////////////////////

  // Initialization of buffer_dump
  Buffers * listBuffers = (Buffers *) malloc(sizeof(Buffers));
  listBuffers->buffer = (BufferOptions **) malloc(sizeof(BufferOptions) * 4);

  listBuffers->buffer[0] = allocateBuffer(bufferSize, bufferFrames,  bufferBytes, "SignalIn");
  listBuffers->buffer[1] = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "SignalOut");
  listBuffers->buffer[2] = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "Autocor");
  listBuffers->buffer[3] = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "FundamentalFrequency");

  // Process tools
  ProcessTools * processTools = (ProcessTools *) malloc(sizeof(ProcessTools));

  processTools->auto_corr = (MY_TYPE *) calloc(bufferFrames, sizeof(MY_TYPE));
  processTools->im_input = (MY_TYPE *) calloc(bufferFrames, sizeof(MY_TYPE));
  processTools->frequencies = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->magnitudes = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->phase = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->autotuneSignal = (MY_TYPE *) calloc(bufferFrames, sizeof(MY_TYPE));
  processTools->oldPhase = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));

  listBuffers->processTools = processTools;

  // Ring Buffer for the autocor function
  listBuffers->ringBufferFreq = createRingBuffer(ringBufferSize);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "[INFO] Start of Recording" << endl;

  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)listBuffers, &options );
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    exit( 1 );
  }

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  try {
    adac.startStream();

    char input;
    std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
    std::cin.get(input);

    // Stop the stream.
    adac.stopStream();
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    goto cleanup;
  }

  cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  /////////////////////////////////////////////// End ///////////////////////////////////////////////

  cout << "[INFO] Writting of records" << endl;

  writeBuffer(listBuffers->buffer[0], PATH_RECORD);
  writeBuffer(listBuffers->buffer[1], PATH_RECORD);
  writeBuffer(listBuffers->buffer[2], PATH_RECORD);
  writeBuffer(listBuffers->buffer[3], PATH_RECORD);

  // Cleaning of dynamic allocations
  deallocateBuffer(listBuffers->buffer[0]);
  deallocateBuffer(listBuffers->buffer[1]);
  deallocateBuffer(listBuffers->buffer[2]);
  deallocateBuffer(listBuffers->buffer[3]);

  free(processTools->auto_corr);
  free(processTools->im_input);
  free(processTools->autotuneSignal);
  free(processTools);

  destroyRingBuffer(listBuffers->ringBufferFreq);

  free(listBuffers);

  cout << "[INFO] End" << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}

/********************************************************
   $Id$
       NAME:
                fft - fast fourier transform
       SYNOPSIS:
                int   fft(x, y, m);

                double   x[];   real part
                double   y[];   imaginary part
                int      m;     data size

                return : success = 0
                         fault   = -1
       Naohiro Isshiki          Dec.1995    modified
********************************************************/

int fft(double *x, double *y, const int m)
{
   int j, lmx, li;
   double *xp, *yp;
   double *sinp, *cosp;
   int lf, lix, tblsize;
   int mv2, mm1;
   double t1, t2;
   double arg;
   int checkm(const int);

   /**************
   * RADIX-2 FFT *
   **************/

   if (checkm(m))
      return (-1);

   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = M_PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }

   lf = maxfftsize / m;
   lmx = m;

   for (;;) {
      lix = lmx;
      lmx /= 2;
      if (lmx <= 1)
         break;
      sinp = _sintbl;
      cosp = _sintbl + maxfftsize / 4;
      for (j = 0; j < lmx; j++) {
         xp = &x[j];
         yp = &y[j];
         for (li = lix; li <= m; li += lix) {
            t1 = *(xp) - *(xp + lmx);
            t2 = *(yp) - *(yp + lmx);
            *(xp) += *(xp + lmx);
            *(yp) += *(yp + lmx);
            *(xp + lmx) = *cosp * t1 + *sinp * t2;
            *(yp + lmx) = *cosp * t2 - *sinp * t1;
            xp += lix;
            yp += lix;
         }
         sinp += lf;
         cosp += lf;
      }
      lf += lf;
   }

   xp = x;
   yp = y;
   for (li = m / 2; li--; xp += 2, yp += 2) {
      t1 = *(xp) - *(xp + 1);
      t2 = *(yp) - *(yp + 1);
      *(xp) += *(xp + 1);
      *(yp) += *(yp + 1);
      *(xp + 1) = t1;
      *(yp + 1) = t2;
   }

   /***************
   * bit reversal *
   ***************/
   j = 0;
   xp = x;
   yp = y;
   mv2 = m / 2;
   mm1 = m - 1;
   for (lmx = 0; lmx < mm1; lmx++) {
      if ((li = lmx - j) < 0) {
         t1 = *(xp);
         t2 = *(yp);
         *(xp) = *(xp + li);
         *(yp) = *(yp + li);
         *(xp + li) = t1;
         *(yp + li) = t2;
      }
      li = mv2;
      while (li <= j) {
         j -= li;
         li /= 2;
      }
      j += li;
      xp = x + j;
      yp = y + j;
   }

   return (0);
}

double *dgetmem(int leng)
{
    return ( (double *)getmem(leng, sizeof(double)) );
}

char *getmem(int leng, unsigned size)
{
    char *p = NULL;

    if ((p = (char *)calloc(leng, size)) == NULL){
        fprintf(stderr, "Memory allocation error !\n");
        exit(3);
    }
    return (p);
}

static int checkm(const int m)
{
   int k;

   for (k = 4; k <= m; k <<= 1) {
      if (k == m)
         return (0);
   }
   fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n",m);

   return (-1);
}
