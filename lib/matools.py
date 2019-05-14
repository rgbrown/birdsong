import wavio # Had to install this using pip - scipy.wavefile cannot handle 24 bit wave files
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

class MicrophoneArray:
    def __init__(self, positions, c_sound=343):
        self.c_sound = c_sound
        self.positions = np.array(positions)
        self.n_microphones = len(positions)

    def generate_synthetic_data(self, x, fs, theta):
        """ Generate an array of synthetic data based on the audio in x. Audio will be shifted according to the array geometry in the microphone array and the speed of sound in there also
        """

        # Compute the requisite delay for each microphone in samples
        delay = self.delay_times(theta)*fs
        
        # Work out the Lagrange interpolation filter details
        N = 6 # lagrange interpolation filter order
        int_delay = np.round(delay - N/2).astype('int')
        frac_delay = delay - int_delay # should be in optimal range

        # make it so the smallest of the integer delays is zero
        int_delay -= min(int_delay)
        max_int_delay = max(int_delay)

        # number of samples in the output
        n_s = len(x) - max_int_delay - N

        # The signal itself will start at max_int_delay
        # Filtering should start at max_int_delay - int_delay
        y = []
        for i in range(self.n_microphones):
            h = lagrange_filter(frac_delay[i], N=N)
            x_slice = x[(max_int_delay - int_delay[i]):]
            y.append(scipy.signal.lfilter(h, 1, x_slice)[N:(N+n_s)])
        return np.array(y)
        
    def load_data(self, base_name):
        data = []
        for i in range(self.n_microphones):
            foo = wavio.read(base_name + '{}.wav'.format(i+1))
            data.append(foo.data.astype('float64').flatten() / (2**(foo.sampwidth*8-1)))
        self.fs = foo.rate
        self.data = np.array(data)

    def get_sample(self, i_start, width):
        return self.data[:,i_start:(i_start+width)]

    def delay_times(self, theta):
        # Compute the delays relative to the origin for a sound coming from direction theta
        d = (np.cos(theta)*self.positions[:,0] + 
             np.sin(theta)*self.positions[:,1])
        return -d/self.c_sound

    def steering_vector(self, theta, f):
        # theta in radians, f in Hz
        tau = self.delay_times(theta)
        a = np.exp(-1j*2*np.pi*f*tau)
        return a

    def get_time(self):
        n_samples = len(self.data[0])
        return np.arange(n_samples)/self.fs
        
    # This might be dodgy
    def get_window(self, t_start, width, theta):
        N = 6 # filter order
        n_s = int(np.round(width*self.fs)) # window size (in samples)
        # d is the number of samples (non-integer) that the signal received at each microphone will be ahead of the
        # signal received at the reference position (origin)
        d = 1/self.c_sound * np.dot(self.positions, np.array([np.cos(theta), np.sin(theta)]))
        t_s_corrected = t_start - d
        y = np.zeros((self.n_microphones, n_s))
        for i in range(self.n_microphones):
            y[i] = shifted_window(self.data[i], self.fx, t_s_corrected[i], width) 
            y[i] -= np.mean(y[i])
            
        t = t_start + (np.arange(n_s)/self.fs)
        return y, t
        
def lagrange_filter(D, N=4):
    #assert((N - 1)/2 <= D and (N + 1)/2 >= D)
    h = np.zeros(N+1)
    for n in range(N+1):
        w = np.setdiff1d(np.arange(N+1), n)
        h[n] = np.prod((D - w) / (n - w))
    return(h)
    #return(np.flip(h, 0))
    #return h

def shifted_window(x, fs, t_start, width, N=6):
    n_s = int(np.round(fs*width))
    i_start = fs*t_start
    int_delay = int(round(i_start - N/2)) # N/2 factor so Lagrange interpolation is in optimal range
    frac_delay = i_start - int_delay
    h = lagrange_filter(frac_delay, N)
    y2 = x[int_delay:(int_delay+n_s + N)]
    # Remove the first N entries - they are junk
    y = scipy.signal.lfilter(h, 1, y2)[N:]
    return y

