import numpy as np
from matplotlib import pyplot as plt
from scipy.io.wavfile import write

# --> Adapt the file name (and path) <--
filename = "./SICOM_BE_Signal_Temps_Reel/autotune_project/output/SignalIn"
#filename = "../../output/SignalIn.bin"

# Read binary file (double format).
x = np.fromfile(filename, dtype=np.double)

# Write signal as a .wav file (useless if signal is not audio).
# The sampling rate is hard coded here, --> change it if necessary <--
write(filename + '.wav', 48000, x)

# Plot the signal. --> Adapt the text on the Figure <--
plt.plot(x)
plt.xlabel('Sample')
plt.ylabel('Input Signal')
plt.title('Evolution of Input Signal')
plt.grid()
plt.show()
