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
N = ma.n_microphones

fs, x = wavfile.read('tui.wav')
theta =  280 * np.pi/180
data = ma.generate_synthetic_data(x, fs, theta)

# Add noise on
sigma = 0.00
for y in data:
    y += np.random.normal(0, sigma, y.shape)

y0 = data[4]

# let's get our head around the FFT - following notation from Lathi Chapter 3 
Ts = 1/fs
T0 = Ts*(len(y0) + 1)
K = len(y0)
freqs = np.arange(K)/T0
X = []
for i in range(ma.n_microphones):
    X.append(np.fft.fft(Ts*data[i]))
X = np.matrix(X)

i0 = 4000
def U_matrix(i):
    Rxx = np.dot(X[:, i], X[:, i].H)
    lam, V = np.linalg.eig(Rxx)
    idx = lam.argsort()[::-1]
    lam = lam[idx]
    V = V[:, idx]
    return V, lam

U0 = U_matrix(i0)[0]

def Ryy(i):
    T = 1/np.sqrt(K) * np.matmul(U0, U_matrix(i)[0].H)
    Y = np.matmul(T, X[:, i])
    return np.matmul(Y, Y.H)

R_coh = np.matrix(np.zeros((N, N))).astype('complex128')
for i in range(K):
    R_coh += Ryy(i)







# there is a peak at index 9937
idx = 9337 
f = freqs[idx]
V, lam = U_matrix(idx)
vv = V[:, 0].flatten()

print('----------')
print('Performing MUSIC at {:.5} Hz'.format(f))
print('-----------------------------')
print('Steering vector subspace check:\n')
print('At the correct angle of {:.3}, '.format(theta*180/np.pi) +
'the real parts of the eigenvalues of R_xx are:')
print('\n'.join('    {:.3}'.format(np.real(l)) for l in lam))
print('\nSteering vector / eigenvector of max eigenvalue:')
print((ma.steering_vector(theta, f) / vv).T)


# Let's do the MUSIC spectrum
En = V[:, 1:] 
n_music = 200
theta_range = np.linspace(0, 2*np.pi, n_music)
P_music = np.zeros(n_music)
for i in range(n_music):
    sv = ma.steering_vector(theta_range[i], f)
    vec = np.dot(En.H, ma.steering_vector(theta_range[i], f))
    P_music[i] = 1/np.linalg.norm(vec)**2

plt.title('$\\theta = {}^\circ$'.format(180*theta/np.pi))
plt.plot(theta_range*180/np.pi, P_music)
plt.show()

# Autofocussing MUSIC time



