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
const int ringBufferSize = 15;
double *_sintbl = 0;
int maxfftsize = 48000;
int samplingFrequency = 512;

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
  assert(sizeAutocor >= 5);

  int max=-1;
  int index=4;

  // TODO: adapt this hard coded variable
  for (int i = 10; i < sizeAutocor - 1; i++){
    if ((autocor[i-1] < autocor[i]) && (autocor[i] > autocor[i+1])){
      if(autocor[i] > max){
        max = autocor[i];
        index = i;
      }
    }
  }

  return index;
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

MY_TYPE process(double * input, int size) {
  // MY_TYPE * auto_corr_value = processTools->auto_corr;
  // MY_TYPE * Im_auto_corr_value = processTools->Im_auto_corr;
  // MY_TYPE * frequencies_value = processTools->frequencies;
  // MY_TYPE * magnitudes_value = processTools->magnitudes;

  MY_TYPE * auto_corr_value = (MY_TYPE *) calloc(size, sizeof(MY_TYPE));
  MY_TYPE * Im_auto_corr_value = (MY_TYPE *) calloc(size, sizeof(MY_TYPE));
  MY_TYPE * frequencies_value = (MY_TYPE *) calloc(size, sizeof(MY_TYPE));
  MY_TYPE * magnitudes_value = (MY_TYPE *) calloc(size, sizeof(MY_TYPE));

  // Search fundament frequency
  demi_auto_corr(input, size, auto_corr_value);
  displayTab("auto corr", auto_corr_value, size);

  int index = extractFundamentalFrequency(auto_corr_value, size);

  if (index == -1){
    fprintf(stderr, "Fundamental frequency not found\n");
    exit(1);
  }

  cout << "index" << index << endl;

  // FFT transformation
  int success = fft(auto_corr_value, Im_auto_corr_value, size);

  if(success == -1) {
    fprintf(stderr, "FFT failed\n");
    exit(1);
  }

  // Setup fundamental frequency
  frequencies_value[0] = (MY_TYPE) (1.0/index) * samplingFrequency;
  magnitudes_value[0] = (MY_TYPE) sqrt((auto_corr_value[index] * auto_corr_value[index]) + (Im_auto_corr_value[index] * Im_auto_corr_value[index]));

  // Magnitude
  int count=2;
  int harmIndex = count * index;

  while (harmIndex < size)
  {
    magnitudes_value[count-1] = sqrt((auto_corr_value[harmIndex] * auto_corr_value[harmIndex])+ (Im_auto_corr_value[harmIndex] * Im_auto_corr_value[harmIndex]));
    frequencies_value[count-1] = (1.0/harmIndex) * samplingFrequency;
    count++;
    harmIndex = count * index;
  }

  displayTab("frequences", frequencies_value, size);
  displayTab("magnitudes", magnitudes_value, size);

  MY_TYPE signal = 0;

  for(int i = 0; i< size; i++){
    signal = signal + 2 * magnitudes_value[i] * cos(2 * M_PI * frequencies_value[i] * i);
  }

  free(auto_corr_value);
  free(Im_auto_corr_value);
  free(frequencies_value);
  free(magnitudes_value);

  return signal;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  Buffers * listBuffers = (Buffers *) data;

  BufferOptions * bufferStructIn = listBuffers[0].buffer;
  BufferOptions * bufferStructOut = listBuffers[1].buffer;

  // ProcessTools * processTools = listBuffers[0].processTools;

  memcpy(outputBuffer, inputBuffer, (size_t) (bufferStructIn->bytes));

  // Record input buffer
  int index = write_buff_dump(
    (double *) inputBuffer,
    bufferStructIn->bufferFrameSize,
    (MY_TYPE *) bufferStructIn->bufferDump,
    bufferStructIn->bufferDumpSize,
    &bufferStructIn->indexBufferDump
  );

  cout << "input" << endl;
  cout << bufferStructIn->bufferDump[0] << " | " <<  bufferStructIn->bufferDump[1] << " | " <<  bufferStructIn->bufferDump[2]<< endl;

  /////////////////// Process ////////////////////////
  MY_TYPE signal = process((double *) inputBuffer, (int) bufferStructIn->bufferFrameSize);

  cout << "signal" << signal << endl;

  MY_TYPE tab[bufferStructIn->bufferFrameSize];

  for (int i=0; i<bufferStructIn->bufferFrameSize; i++){
    if(signal > 1){
      tab[i] = 1;
    }
    else if(signal < -1){
      tab[i] = -1;
    }
    else
      tab[i] = signal;
  }

  ////////////////////////////////////////////////////

  // Record output buffer
  int index2 = write_buff_dump(
    (MY_TYPE *) tab,
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
  Buffers * listBuffers = (Buffers *) malloc(sizeof(listBuffers) * 3);
  BufferOptions * bufferIn = allocateBuffer(bufferSize, bufferFrames,  bufferBytes, "SignalIn");
  BufferOptions * bufferOut = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "SignalOut");
  //BufferOptions * fundamentalFrequency = allocateBuffer(bufferSize, bufferFrames , bufferBytes, "FundamentalFrequency");

  listBuffers[0].buffer = bufferIn;
  listBuffers[1].buffer = bufferOut;
  //listBuffers[2].buffer = fundamentalFrequency;

  ProcessTools * processTools = (ProcessTools *) malloc(sizeof(ProcessTools));

  processTools->auto_corr = (MY_TYPE *) malloc(sizeof(MY_TYPE) * bufferFrames);
  processTools->Im_auto_corr = (MY_TYPE *) malloc(sizeof(MY_TYPE) * bufferFrames);
  processTools->magnitudes = (MY_TYPE *) malloc(sizeof(MY_TYPE) * bufferFrames);
  processTools->frequencies = (MY_TYPE *) malloc(sizeof(MY_TYPE) * bufferFrames);

  listBuffers[0].processTools = processTools;
  listBuffers[1].processTools = processTools;
  //listBuffers[2].processTools = processTools;

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
  //writeBuffer(fundamentalFrequency, PATH_RECORD);

  // Cleaning of dynamic allocations
  deallocateBuffer(bufferIn);
  deallocateBuffer(bufferOut);
  //deallocateBuffer(fundamentalFrequency);

  free(listBuffers);
  // free(processTools->auto_corr);
  // free(processTools->Im_auto_corr);
  // free(processTools->magnitudes);
  // free(processTools->frequencies);
  // free(processTools);

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
