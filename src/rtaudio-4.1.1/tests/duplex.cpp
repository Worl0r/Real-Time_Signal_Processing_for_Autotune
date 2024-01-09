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
#include "somefunc.cpp"

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

/////////////////////////////////////////////// Configuration ///////////////////////////////////////////////

const string PATH_RECORD = "../../../files/";
int AUTOTUNE = false;
int DISPLAY = false;
unsigned int bufferFrames = 512;
const int bufferSize = 500000;
const int ringBufferSize = 2;
double samplingFrequency = 48000;
int nbrHarmonics = bufferFrames / 2;
int nbrDemiTons = 12;
int nfftSize = bufferFrames;
int jumpedIdx = 4;

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
    // We skip zero values
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
    fprintf(stderr, "[ERROR] Memory allocation error for BufferOptions\n");
    exit(1);
  }

  bufferOptions->bufferDumpSize = bufferDumpSize;

  bufferOptions->bufferDump = (MY_TYPE*) calloc(bufferOptions->bufferDumpSize, sizeof(MY_TYPE));
  if (bufferOptions->bufferDump == NULL) {
    fprintf(stderr, "[ERROR] Memory allocation error for bufferDump\n");
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

/////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////// process //////////////////////////////////////////////

// This function extract the funcdamental frequency from the autocor
int extractFundamentalFrequency(MY_TYPE * autocor, int sizeAutocor){

  assert(autocor);
  assert(sizeAutocor >= 5);

  MY_TYPE max = 0.0;
  int idx = 0;

  for(int i = jumpedIdx; i < (sizeAutocor-1);i ++){
    if(autocor[i-1] < autocor[i] && autocor[i] > autocor[i+1] && max < autocor[i]){
      max = autocor[i];
      idx = i;
    }
  }

  return idx;
}

// Compute the autocorrelation of the input signal
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

// Change the fundamental frequency to the nearest frequency in nbrDemiTons
MY_TYPE changeFundamentalFrequency(MY_TYPE fundamentalFrequency, int nbrDemiTons, int size){
  fundamentalFrequency = 12.0 * log2(fundamentalFrequency);

  fundamentalFrequency = round(fundamentalFrequency);
  fundamentalFrequency =((int) fundamentalFrequency / (int) nbrDemiTons) * nbrDemiTons;

  fundamentalFrequency = pow(2, fundamentalFrequency / 12.0);

  return fundamentalFrequency;
}

// Compute the fundamental frequency and save it in the ring buffer
int saveFundamentalFrequency(RingBuffer * ringBuffer, ProcessTools * processTools, double * input, int size){

  // Clear autocorrelation table
  memset(processTools->autoCorr, 0, size * sizeof(MY_TYPE));

  // Compute autocorrelation
  demi_auto_corr(input, size, processTools->autoCorr);

  int index = extractFundamentalFrequency(processTools->autoCorr, size);

  // Save fundamental frequency
  writeRingBuffer(ringBuffer, index);

  return 1;
}

// This function apply the autotune process
void autotuneProcess(double * input, RingBuffer * ringBufferFreq, ProcessTools * processTools, int size) {

  // Compute fundamental frequency
  int val = saveFundamentalFrequency(ringBufferFreq, processTools, input, size);

  if (val == -1){
    fprintf(stderr, "[ERROR] Fundamental frequency not found\n");
    exit(1);
  }

  MY_TYPE im_input[size];
  for (int i = 0; i < size; i++){
    im_input[i] = 0;
  }

  // FFT transformation
  int success = fftr(input, im_input, nfftSize);

  if(success == -1) {
    fprintf(stderr, "[ERROR] FFT failed\n");
    exit(1);
  }

  // Compute mean of fundamental frequency through the ring buffer
  int meanFreq = meanRingBuffer(ringBufferFreq);

  if(meanFreq == 0){
    cout << "[INFO] Fundamental frequency not found" << endl;
    return;
  }

  // Convert meanFreq to Hz
  MY_TYPE meanFreqHz = (MY_TYPE) (samplingFrequency / meanFreq);

  // Transform fundamental frequency
  if (AUTOTUNE){
    meanFreqHz = changeFundamentalFrequency(meanFreqHz, nbrDemiTons, nbrHarmonics);
  }

  // Adapt the fundamental frequency to the size of the buffer
  meanFreq = (int)round(meanFreqHz * size / (samplingFrequency));

  // Setup fundamental frequency and magnitude
  processTools->frequencies[0] = meanFreqHz;
  processTools->magnitudes[0] = (MY_TYPE) sqrt((input[meanFreq] * input[meanFreq]) + (im_input[meanFreq] * im_input[meanFreq]));
  processTools->magnitudes[0] = processTools->magnitudes[0] / nfftSize;
  processTools->phase[0] = fmod(processTools->oldPhase[0] + 2.0 *  M_PI * processTools->oldFrequencies[0] * ((double) size / samplingFrequency), 2.0 * M_PI);

  // Idem for other harmonics
  for (int j = 2; j * meanFreq < size /2; j++){
    processTools->frequencies[j-1] = j * meanFreqHz;
    processTools->magnitudes[j-1] = sqrt((input[j * meanFreq] * input[j * meanFreq]) + (im_input[j * meanFreq] * im_input[j * meanFreq]));
    processTools->magnitudes[j-1] = processTools->magnitudes[j-1] / nfftSize;
    processTools->phase[j-1] = fmod(processTools->oldPhase[j-1] + 2.0 *  M_PI * processTools->oldFrequencies[j-1] * ((double) size / samplingFrequency), 2.0 * M_PI);
  }

  // Update old phase and old frequencies
  memcpy(processTools->oldPhase, processTools->phase, nbrHarmonics * sizeof(MY_TYPE));
  memcpy(processTools->oldFrequencies, processTools->frequencies, nbrHarmonics * sizeof(MY_TYPE));

  // Compute the transformed signal through the time
  for(int i = 0; i< size; i++){
    MY_TYPE sum = 0.0;

    // Synthesis of the signal
    for(int j = 0; j * meanFreq < size / 2; j++){
      sum = sum + 2.0 * processTools->magnitudes[j] * cos(2.0 * M_PI * (double) (i+1) * (meanFreqHz * (j+1) / samplingFrequency) + processTools->phase[j]);
    }

    // Thresholding
    if(sum > 1){
      processTools->autotuneSignal[i] = 1;
    }
    else if(sum < -1){
      processTools->autotuneSignal[i] = -1;
    }
    else{
      processTools->autotuneSignal[i] = sum;
    }
  }

  // Display
  displayRingBuffer(ringBufferFreq);

  if (DISPLAY){
    cout << "[DEBUG] frequencies_value[0] in Hz after transformation: " << meanFreqHz << endl;
    cout << "[DEBUG] frequencies_value[0] in idx after transformation: " << meanFreq << endl;
    displayTab("Magnitudes", processTools->magnitudes, nbrHarmonics);
    displayTab("Frequencies", processTools->frequencies, nbrHarmonics);
    // displayTab("oldPhase", processTools->oldPhase, nbrHarmonics);
    // displayTab("Phase", processTools->phase, nbrHarmonics);
    //displayTab("SignalTab", (processTools->autotuneSignal), size);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  // Variable extraction
  Buffers * listBuffers = (Buffers *) data;

  BufferOptions * bufferStructIn = listBuffers->buffer[0];
  BufferOptions * bufferStructOut = listBuffers->buffer[1];
  BufferOptions * autocor = listBuffers->buffer[2];
  BufferOptions * fundamentalFrequency = listBuffers->buffer[3];

  RingBuffer * ringBufferFreq = listBuffers->ringBufferFreq;

  ProcessTools * processTools = listBuffers->processTools;

  // Record input buffer for figures
  int index = write_buff_dump(
    (double *) inputBuffer,
    bufferStructIn->bufferFrameSize,
    (MY_TYPE *) bufferStructIn->bufferDump,
    bufferStructIn->bufferDumpSize,
    &bufferStructIn->indexBufferDump
  );

  /////////////////// Autotune Process //////////////////////

  autotuneProcess((double *) inputBuffer, ringBufferFreq, processTools, (int) bufferStructIn->bufferFrameSize);

  ///////////////////////////////////////////////////////////

  // Copy input buffer to output buffer to hear the effects of the processing in realtime
  memcpy(outputBuffer, (void *) processTools->autotuneSignal, (size_t) (bufferStructIn->bytes));

  // Record output buffer for figures
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

  cout << "[INFO] Dynamic memory is being allocated..." << endl;

  // Initialization of buffer_dump
  Buffers * listBuffers = (Buffers *) malloc(sizeof(Buffers));
  listBuffers->buffer = (BufferOptions **) calloc(4, sizeof(BufferOptions*));

  listBuffers->buffer[0] = allocateBuffer(bufferSize, bufferFrames,  bufferBytes, "SignalIn");
  listBuffers->buffer[1] = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "SignalOut");
  listBuffers->buffer[2] = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "Autocor");
  listBuffers->buffer[3] = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "FundamentalFrequency");

  // Process tools
  ProcessTools * processTools = (ProcessTools *) malloc(sizeof(ProcessTools));

  processTools->autoCorr = (MY_TYPE *) calloc(bufferFrames, sizeof(MY_TYPE));
  processTools->frequencies = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->oldFrequencies = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->magnitudes = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->phase = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));
  processTools->autotuneSignal = (MY_TYPE *) calloc(bufferFrames, sizeof(MY_TYPE));
  processTools->oldPhase = (MY_TYPE *) calloc(nbrHarmonics, sizeof(MY_TYPE));

  listBuffers->processTools = processTools;

  // Ring Buffer for the fundamental frequency
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

  cout << "[INFO] Writting of records..." << endl;

  writeBuffer(listBuffers->buffer[0], PATH_RECORD);
  writeBuffer(listBuffers->buffer[1], PATH_RECORD);
  writeBuffer(listBuffers->buffer[2], PATH_RECORD);
  writeBuffer(listBuffers->buffer[3], PATH_RECORD);

  // Cleaning of dynamic allocations
  cout << "[INFO] Cleaning of the dynamic memory..." << endl;
  deallocateBuffer(listBuffers->buffer[0]);
  deallocateBuffer(listBuffers->buffer[1]);
  deallocateBuffer(listBuffers->buffer[2]);
  deallocateBuffer(listBuffers->buffer[3]);

  free(processTools->autoCorr);
  free(processTools->autotuneSignal);
  free(processTools->frequencies);
  free(processTools->oldFrequencies);
  free(processTools->magnitudes);
  free(processTools->phase);
  free(processTools->oldPhase);
  free(processTools);

  destroyRingBuffer(listBuffers->ringBufferFreq);

  free(listBuffers->buffer);
  free(listBuffers);

  cout << "[INFO] End" << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}

