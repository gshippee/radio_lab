import ugradio
import matplotlib.pylab as plt
import numpy as np
import argparse
import os

sample_rate = 1024
dt = 1.0/sample_rate
t = np.arange(sample_rate)*dt  # 1 second of samples
freq = 5
amp = 1.0
sine = amp*np.sin(2*np.pi*freq*t)

def hamming_window(data, num_samples):
	for n in range(num_samples):
		current_window_value = 0.54-0.46*np.cos(2*np.pi*n/(num_samples-1))
		data[n] = current_window_value*data[n]
	return data

ones = np.ones(1024)
ones_hamming = hamming_window(ones, 1024)

plt.plot(ones_hamming)

def plot_freq_data(power, sample_rate, xlabel = 'Frequency (Hz)', ylabel = 'Power (V$\cdot$s)'):
	FT = np.abs(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(power))))**2
	freqs = np.fft.fftshift(np.fft.fftfreq(len(power), 1./sample_rate))
	freq_res = round(sample_rate/len(power)*1000,2)
	print('Max power found at: ' + str(freqs[np.argmax(FT)]) + ' Hz')
	plt.figure()
	plt.plot(freqs, FT)
	plt.xlabel(xlabel, fontsize = 14)
	plt.ylabel(ylabel, fontsize = 14)
	plt.xlim(0,10)
	plt.title('Power Spectrum $\Delta$f = ' + str(freq_res) + ' kHz', fontsize = 18) 


plt.show()
