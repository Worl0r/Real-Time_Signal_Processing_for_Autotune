import numpy as np
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
import os
import math

## Options
filename = "./SICOM_BE_Signal_Temps_Reel/autotune_project/files/"
sampling_rate = 48000


def plotSignal(name, frequence, color):
    # Read binary file (double format).
    y = np.fromfile(os.path.join(filename, name), dtype=np.double)

    # Write signal as a .wav file (useless if signal is not audio).
    # The sampling rate is hard coded here, --> change it if necessary <--
    write(os.path.join(filename, name) + '.wav', frequence, y)

    # Plot the signal.
    time = np.arange(y.shape[0]) / frequence  # Create a time array
    plt.figure()
    plt.plot(time, y, color)

    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Evolution of ' + name + ' Over Time')
    plt.grid()
    plt.savefig(os.path.join(filename, name))

def plotCompareSignal(name1, name2, frequence, color1, color2):
    # Read binary file (double format).
    y1 = np.fromfile(os.path.join(filename, name1), dtype=np.double)
    y2 = np.fromfile(os.path.join(filename, name2), dtype=np.double)

    # Plot the signal.
    time = np.arange(y1.shape[0]) / frequence  # Create a time array
    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    ax[0].plot(time, y1, color1)
    ax[1].plot(time, y2, color2)

    fig.suptitle('Evolution of ' + name1 + ' and ' + name2 + ' Over Time')
    #define the axes
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_title(name1)
    ax[0].grid()

    if name2 != "FundamentalFrequency":
        ax[1].set_ylabel('Amplitude')
    else:
        mask = np.isfinite(y2)
        y2 = y2[mask]
        ax[1].set_ylabel('Frequency (Hz)')
        if np.max(y2) > 700:
            ax[1].set_yticks(np.arange(0, np.max(y2), 500))
        else:
            ax[1].set_yticks(np.arange(0, np.max(y2), 50))


    ax[1].set_xlabel('Time (s)')

    ax[1].set_title(name2)
    ax[1].grid()

    plt.savefig(os.path.join(filename, name1 + name2))

def computeFundamentalFrequency2(name, frequence):
    y = np.fromfile(os.path.join(filename, name), dtype=np.double)
    fft_spectrum = np.fft.rfft(y)
    return np.fft.rfftfreq(y.size, d=1./frequence)[np.argmax(np.abs(fft_spectrum))]

def main():
    plotSignal("SignalIn", sampling_rate, 'b')
    plotSignal("SignalOut",sampling_rate, 'r')
    plotSignal("Autocor", sampling_rate, 'r')
    plotSignal("FundamentalFrequency", sampling_rate, 'r')

    f = computeFundamentalFrequency2("SignalOut", sampling_rate)
    print("Fundamental Frequency: " + str(f))

    plotCompareSignal("SignalIn", "SignalOut", sampling_rate, 'b', 'r')
    plotCompareSignal("SignalIn", "FundamentalFrequency", sampling_rate, 'b', 'r')


    print("[END] Finished")


if __name__ == '__main__':
    main()
