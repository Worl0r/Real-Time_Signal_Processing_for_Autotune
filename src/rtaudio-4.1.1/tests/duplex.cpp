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

RingBuffer::RingBuffer(int size): sizeRingBuffer(size), indexRead(-1), indexWrite(0){
  dataRingBuffer = new int[sizeRingBuffer];
};

RingBuffer::~RingBuffer(){
  free(dataRingBuffer);
};

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

const string PATH_RECORD = "./SICOM_BE_Signal_Temps_Reel/autotune_project/output/";
const int bufferSize = 10000;

/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////// Utils ///////////////////////////////////////////////

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
}

BufferOptions * allocateBuffer(int bufferSize, unsigned int bufferBytes, int indexBufferDump = 0){
  BufferOptions * bufferOptions = (BufferOptions *) malloc(sizeof(BufferOptions));
  if (bufferOptions == NULL) {
    fprintf(stderr, "Memory allocation error for BufferOptions\n");
    exit(1);
  }

  bufferOptions->bufferSize = 100000;
  bufferOptions->bufferDump = (MY_TYPE*) calloc(bufferOptions->bufferSize, sizeof(MY_TYPE));
  if (bufferOptions->bufferDump == NULL) {
    fprintf(stderr, "Memory allocation error for bufferDump\n");
    exit(1);
  }

  bufferOptions->indexBufferDump = 0;
  bufferOptions->bytes = bufferBytes;

  return bufferOptions;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data )
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  //unsigned int *bytes = (unsigned int *) data;
  BufferOptions *bufferOptions = (BufferOptions *) data;

  memcpy( outputBuffer, inputBuffer, (size_t) (bufferOptions->bytes));

  int index = write_buff_dump(
    (double *) outputBuffer,
    bufferOptions->bytes,
    (MY_TYPE *) bufferOptions->bufferDump,
    bufferOptions->bufferSize,
    &bufferOptions->indexBufferDump
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

  BufferOptions * bufferIn = allocateBuffer(bufferSize, bufferBytes);
  BufferOptions * bufferOut = allocateBuffer(bufferSize, bufferBytes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)bufferIn, &options );
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

deallocateBuffer(bufferIn);

////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}


