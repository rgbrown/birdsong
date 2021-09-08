import matools
import numpy as np
import matplotlib.pyplot as plt

# Set up the problem data
positions = 1e-3*np.array([
    [0, 0],
    [42.6, -0.5], 
    [21.5, -37.7],
    [-21.2, -37.7], 
    [-42.3, -0.8],
    [-21.2, 35.8],
    [21.4, 36.2],
    ])

fs = 44100
ma = matools.MicrophoneArray(positions)
basename = "../datasets/sample_alberto_uma8/Audio Track-"
ma.load_data_wavfile(basename)

# Set up a window
i_start = 8132000
i_finish = 8148000

data_slice = ma.data[:,i_start:i_finish]
