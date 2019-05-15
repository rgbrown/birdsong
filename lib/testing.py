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
theta =  90 * np.pi/180
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
idx = 9937
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

# Let's do the MUSIC spectrum
En = V[:, np.setdiff1d(np.arange(5), idx)]
n_music = 200
theta_range = np.linspace(0, 2*np.pi, n_music)
P_music = np.zeros(n_music)
for i in range(n_music):
    sv = ma.steering_vector(theta_range[i], f)
    vec = np.dot(En.H, ma.steering_vector(theta_range[i], f))
    P_music[i] = 1/np.linalg.norm(vec)**2

plt.plot(theta_range*180/np.pi, P_music)
plt.show()

