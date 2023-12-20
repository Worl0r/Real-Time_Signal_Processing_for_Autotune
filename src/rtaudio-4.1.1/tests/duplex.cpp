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

/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////// Utils ///////////////////////////////////////////////

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

void derivative(MY_TYPE * f, int size, MY_TYPE *derivative){
  derivative[0] = 0;
  for (int i = 1; i < size-1; i++){
    derivative[i] = (f[i + 1] - f[i - 1]) / 2;
  }
}

// The goal is to extract the second maximum for the autocor function.
int extractFundamentalFrequency(MY_TYPE * autocor, int sizeAutocor){
  assert(autocor);

  // // First derivative
  // MY_TYPE * firstDerivative = (MY_TYPE *) malloc(sizeof(MY_TYPE) * (sizeAutocor-1));
  // derivative(autocor, sizeAutocor, firstDerivative);

  // // Second derivative
  // MY_TYPE * secondDerivative = (MY_TYPE *) malloc(sizeof(MY_TYPE) * (sizeAutocor-2));
  // derivative(firstDerivative, sizeAutocor, secondDerivative);

  // // Keep the first maximum because we cannot compute a zero derivative for the first index.
  // for (int i = 1; i < sizeAutocor-2; i++){
  //   if ((firstDerivative[i] >= 0 && firstDerivative[i+1] <= 0) ||
  //     (firstDerivative[i] <= 0 && firstDerivative[i+1] >= 0)
  //   ){
  //     if (secondDerivative[i] <= 0){
  //       free(firstDerivative);
  //       free(secondDerivative);
  //       // [IMPORTANT] The index stats to zero.
  //       return i;
  //     }
  //   }
  // }

  int * maxList = (int *) malloc(sizeof(int) * (sizeAutocor));
  int count = 0;

  for (int i = 1; i < sizeAutocor - 1; i++){
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

/////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////// process //////////////////////////////////////////////

void process(double * input, int size, BufferOptions * fundamentalFrequency) {

  MY_TYPE* auto_corr = (MY_TYPE*) calloc(size ,sizeof(MY_TYPE));

  demi_auto_corr(input, size, auto_corr);

  int index = extractFundamentalFrequency(auto_corr, size);

  free(auto_corr);

  MY_TYPE tab[fundamentalFrequency->bufferFrameSize];

  for (int i=0; i<fundamentalFrequency->bufferFrameSize; i++){
    tab[i] = (MY_TYPE) index;
  }

  int _ = write_buff_dump(
    tab,
    fundamentalFrequency->bufferFrameSize,
    fundamentalFrequency->bufferDump,
    fundamentalFrequency->bufferDumpSize,
    &fundamentalFrequency->indexBufferDump
  );

}
/////////////////////////////////////////////////////////////////////////////////////////////////////

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  // Unsigned int *bytes = (unsigned int *) data;
  Buffers *listBuffers = (Buffers *) data;
  BufferOptions * bufferStructIn = listBuffers[0].buffer;
  BufferOptions * bufferStructOut = listBuffers[1].buffer;
  BufferOptions * fundamentalFrequency = listBuffers[2].buffer;

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

  process((double *) inputBuffer, (int) bufferStructIn->bytes, fundamentalFrequency);

  ////////////////////////////////////////////////////

  // Record output buffer
  int index2 = write_buff_dump(
    (double *) outputBuffer,
    bufferStructOut->bufferFrameSize,
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

  Buffers * listBuffers = (Buffers *) malloc(sizeof(listBuffers) * 3);
  BufferOptions * bufferIn = allocateBuffer(bufferSize, bufferFrames,  bufferBytes, "SignalIn");
  BufferOptions * bufferOut = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "SignalOut");
  BufferOptions * fundamentalFrequency = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "FundamentalFrequency");
  listBuffers[0].buffer = bufferIn;
  listBuffers[1].buffer = bufferOut;
  listBuffers[2].buffer = fundamentalFrequency;

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
writeBuffer(bufferIn, PATH_RECORD);
writeBuffer(bufferOut, PATH_RECORD);
writeBuffer(fundamentalFrequency, PATH_RECORD);

// Cleaning of dynamic allocations
deallocateBuffer(bufferIn);
deallocateBuffer(bufferOut);
deallocateBuffer(fundamentalFrequency);
free(listBuffers);

// Tests
cout << "[INFO] [TEST] Test of Ring Buffer" << endl;
RingBuffer ringBuffer = RingBuffer(4);
ringBuffer.displayRingBuffer();
ringBuffer.writeRingBuffer((MY_TYPE) 30);
ringBuffer.displayRingBuffer();
ringBuffer.writeRingBuffer((MY_TYPE) 40);
ringBuffer.writeRingBuffer((MY_TYPE) 50);
ringBuffer.displayRingBuffer();

MY_TYPE value = ringBuffer.readRingBuffer();
cout << "My value is :" << value << endl;
ringBuffer.displayRingBuffer();
ringBuffer.writeRingBuffer((MY_TYPE) 90);
ringBuffer.displayRingBuffer();
ringBuffer.writeRingBuffer((MY_TYPE) 100);
ringBuffer.displayRingBuffer();

cout << "\n" << endl;

cout << "[INFO] [TEST] Test Maximum" << endl;

int length = 100;
int * tmp = (int *) malloc(sizeof(int) * length);
MY_TYPE * f = (MY_TYPE *) malloc(sizeof(MY_TYPE) * length);
float fre = 1.0/5.0;

for (int i=0; i<length; i++){
  tmp[i] = i;
}
for (int i=0; i<length; i++){
  f[i] = cos(tmp[i] * 2 * M_PI * fre);
}

int i = extractFundamentalFrequency(f, length);
cout << "\nThe period is : " << i+1 << " indexs." << endl;

cout << "\n" << endl;

cout << "[INFO] End" << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}

////////////////////////////////////////////// Ring Buffer /////////////////////////////////////////

// Constructor and Desstructor
RingBuffer::RingBuffer(int size): sizeRingBuffer(size), indexRead(-1), indexWrite(0){
  buffer = new MY_TYPE[sizeRingBuffer];

  for (int i = 0; i < sizeRingBuffer; i++){
    buffer[i] = 0;
  }
}

RingBuffer::~RingBuffer(){
  free(buffer);
}

// Methods
int RingBuffer::writeRingBuffer(int data){
  if (data == NULL){

    cout << "[ERROR] Data is NULL in writeRingBuffer" << endl;
    return -1;
  }

  if (indexWrite == indexRead ){
    indexRead++;
    indexRead = indexRead % sizeRingBuffer;
  }

  if (indexRead == -1){
    indexRead = 0;
  }

  buffer[indexWrite] = data;
  indexWrite++;
  indexWrite = indexWrite % sizeRingBuffer;

  return 0;
}

MY_TYPE RingBuffer::readRingBuffer(){

    if (indexRead == -1){
      cout << "[ERROR] There is nothing to read in the ring buffer in readRingBuffer" << endl;
      return NULL;
    }

    MY_TYPE data = buffer[indexRead];
    indexRead++;
    indexRead = indexRead % sizeRingBuffer;
    return data;
}

void RingBuffer::displayRingBuffer(){

  std::time_t result = std::time(nullptr);
  cout << "[LOG][TEST] You display your ring buffer at " << std::asctime(std::localtime(&result)) << endl;
  cout << "[ ";

  for(int i=0; i < sizeRingBuffer; i++){
    if (i == sizeRingBuffer-1){
      cout << buffer[i];
    }
    else{
      cout << buffer[i] << " | ";
    }
  }

  cout << " ]" << endl;
  cout << "Read index is : " << indexRead << endl;
  cout << "Write index is : " << indexWrite << endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
