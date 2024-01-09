import numpy as np
from matplotlib import pyplot as plt
from scipy.io.wavfile import write
import os
import math

# --> Adapt the file name (and path) <--
filename = "./SICOM_BE_Signal_Temps_Reel/autotune_project/files/"

def plotSignal(name, frequence, color):
    # Read binary file (double format).
    x = np.fromfile(os.path.join(filename, name), dtype=np.double)

    # Write signal as a .wav file (useless if signal is not audio).
    # The sampling rate is hard coded here, --> change it if necessary <--
    write(os.path.join(filename, name) + '.wav', frequence, x)

    # Plot the signal. --> Adapt the text on the Figure <--
    plt.figure()
    if name == "FundamentalFrequency":
        plt.plot(((1 / x) * frequence).astype(float), color)
    else:
        plt.plot(x, color)

    plt.xlabel('Sample')
    plt.ylabel(name)
    plt.title('Evolution of ' + name)
    plt.grid()
    plt.savefig(os.path.join(filename, name))

def main():
    plotSignal("SignalIn", 48000, 'b')
    plotSignal("SignalOut", 48000, 'r')
    plotSignal("Autocor", 48000, 'r')
    plotSignal("FundamentalFrequency", 48000, 'r')


    print("[END] Finished")


if __name__ == '__main__':
    main()
