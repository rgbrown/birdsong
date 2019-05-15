import matools
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

# Set up problem data
dx = 0.05 
ma = matools.MicrophoneArray([
    [dx, 0], 
    [0, dx], 
    [-dx, 0], 
    [0, -dx], 
    [0,0]])
N = ma.n_microphones
fs, x = wavfile.read('tui.wav')
theta = 280 * np.pi/180
data = ma.generate_synthetic_data(x, fs, theta)

# Add noise on
sigma = 0.005
for y in data:
    y += np.random.normal(0, sigma, y.shape)

y0 = data[4]

# let's get our head around the FFT - following notation from Lathi Chapter 3 
Ts = 1/fs
K = len(y0)
T0 = Ts*(K + 1) # T0 is the end time
freqs = np.arange(K)/T0 
X = [] # Matrix of Fourier transforms, one per row
for i in range(ma.n_microphones):
    X.append(np.fft.fft(Ts*data[i]))
X = np.matrix(X)

i0 = 9337 # Index of reference frequency for AF-MUSIC
def eig_sorted(A):
    lam, V = np.linalg.eig(A)
    idx = lam.argsort()[::-1]
    lam = lam[idx]
    V = V[:, idx]
    return lam, V

def U_matrix(i):
    Rxx = np.dot(X[:, i], X[:, i].H)
    lam, V = eig_sorted(Rxx)
    return V



def music(idx, n_music=200):
    """ Perform MUSIC at frequency freqs[idx]"""
    f = freqs[idx]
    Rxx = np.dot(X[:, idx], X[:, idx].H)
    lam, V = eig_sorted(Rxx)
    En = V[:, 1:] # Noise subspace for one source

    theta_range = np.linspace(0, 2*np.pi, n_music)
    P_music = np.zeros(n_music)
    for i in range(n_music):
        sv = ma.steering_vector(theta_range[i], f)
        vec = np.dot(En.H, ma.steering_vector(theta_range[i], f))
        P_music[i] = 1/np.linalg.norm(vec)**2

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
    return P_music, theta_range

def afmusic(i0, n_music=200):
    f0 = freqs[i0]
    U0 = U_matrix(i0)
    def Ryy(i):
        T = 1/np.sqrt(K) * np.matmul(U0, U_matrix(i).H)
        Y = np.matmul(T, X[:, i])
        return np.matmul(Y, Y.H)
         
    R_coh = np.matrix(np.zeros((N, N))).astype('complex128')
    for i in range(K):
        R_coh += Ryy(i)

    lam, V = eig_sorted(R_coh)
    En = V[:, 1:] 

    theta_range = np.linspace(0, 2*np.pi, n_music)
    P_afmusic = np.zeros(n_music)
    for i in range(n_music):
        sv = ma.steering_vector(theta_range[i], f0)
        vec = np.dot(En.H, ma.steering_vector(theta_range[i], f0))
        P_afmusic[i] = 1/np.linalg.norm(vec)**2
    return P_afmusic, theta_range

# Let's do the MUSIC spectrum
i_music = 9337
P_music, theta_range = music(i_music)
plt.figure(1)
plt.title('MUSIC at {:.4} Hz: $\\theta_0 = {}^\circ$'.format(freqs[i_music], 180*theta/np.pi))
plt.plot(theta_range*180/np.pi, P_music)
plt.xlabel('$\\theta$')

# Autofocussing MUSIC time
i0 = 4000
P_afmusic, theta_range = afmusic(i0)
plt.figure(2)
plt.title('AF_MUSIC with $f_0 = {:.5}$ Hz: $\\theta_0 = {}^\circ$'.format(freqs[i0], 180*theta/np.pi))
plt.plot(theta_range*180/np.pi, P_afmusic)
plt.xlabel('$\\theta$')
plt.show()


