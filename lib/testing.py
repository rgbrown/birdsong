import matools
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

ma = matools.MicrophoneArray([
    [0.01, 0], 
    [0, 0.01], 
    [-0.01, 0], 
    [0, -0.01], 
    [0,0]])

fs, x = wavfile.read('tui.wav')
theta = np.pi/3
y = ma.generate_synthetic_data(x, fs, theta)

y0 = y[4]

# let's get our head around the FFT - following notation from Lathi Chapter 3
Ts = 1/fs
T0 = Ts*(len(y0) + 1)
N0 = len(y0)
freqs = np.arange(N0)/T0
Y = []
for i in range(ma.n_microphones):
    Y.append(np.fft.fft(Ts*y[i]))
Y = np.matrix(Y)

# there is a peak at index 9937
idx = 4000
f = freqs[idx]
Rxx = np.dot(Y[:,idx], Y[:,idx].H)
lam, V = np.linalg.eig(Rxx)
idx = np.argmax(abs(lam))
vv = V[:, idx].flatten()

print('----------')
print('At frequency {:.5} Hz, there is a peak'.format(f))
print('The absolute values of the eigenvalues of R_xx are:')
print(np.abs(lam))
print('Steering vector / eigenvector of max eigenvalue:')
print((ma.steering_vector(theta, f) / vv).T)

