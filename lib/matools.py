import wavio # Had to install this using pip - scipy.wavefile cannot handle 24 bit wave files
import numpy as np
import scipy.signal
class MicrophoneArray:
    def __init__(self, positions, c_sound=343):
        self.c_sound = c_sound
        self.positions = np.array(positions)
        self.n_microphones = len(positions)
        
    def load_data(self, base_name):
        data = []
        for i in range(self.n_microphones):
            foo = wavio.read(base_name + '{}.wav'.format(i+1))
            data.append(foo.data.astype('float64').flatten() / (2**(foo.sampwidth*8-1)))
        self.fs = foo.rate
        self.data = np.array(data)

    def get_time(self):
        n_samples = len(self.data[0])
        return np.arange(n_samples)/self.fs
        
    def get_window(self, t_start, width, theta):
        N = 6 # filter order
        n_s = int(np.round(width*self.fs)) # window size (in samples)
        # d is the number of samples (non-integer) that the signal received at each microphone will be ahead of the
        # signal received at the reference position (origin)
        d = self.fs/self.c_sound * np.dot(self.positions, np.array([np.cos(theta), np.sin(theta)]))
        t_s_corrected = self.fs*t_start - d
        y = np.zeros((self.n_microphones, n_s))
        
        for i in range(self.n_microphones):
            t_s = t_s_corrected[i]
            int_delay = int(round(t_s - N/2)) # the N/2 factor is so that the interpolation is performed in the optimal range
            frac_delay = t_s - int_delay
            h = lagrange_filter(frac_delay, N)
            y2 = self.data[i, int_delay:(int_delay + n_s + N)]
            # first N are junk entries, so strip em off
            y[i] = scipy.signal.lfilter(h, 1, y2)[N:] 
            y[i] -= np.mean(y[i])
        t = t_start + (np.arange(n_s)/self.fs)
        return y, t
        
def lagrange_filter(D, N=4):
    #assert((N - 1)/2 <= D and (N + 1)/2 >= D)
    h = np.zeros(N+1)
    for n in range(N+1):
        w = np.setdiff1d(np.arange(N+1), n)
        h[n] = np.prod((D - w) / (n - w))
    return(np.flip(h, 0))
    #return h
