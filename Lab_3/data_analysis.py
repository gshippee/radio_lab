import ugradio
import matplotlib.pylab as plt
import numpy as np
import argparse
import os

sun = False
crab = True

if crab == True:
	d = np.load('CrabDataR40.npz')
	ra = 1.45967339498
	dec = 0.38422550818
	chunk_divisor = 186
#	chunk_divisor = 124

if sun == True:
	d = np.load('SunDataR51.npz')
	ra = 0
	chunk_divisor = 125

time = d['arr_0']
print(len(time))
data = d['arr_1']

time = time[0:chunk_divisor**2]
data = data[0:chunk_divisor**2]


def remove_DC(power, chunk_divisor):
	power = power[:int(len(power)/chunk_divisor)*chunk_divisor]
	chunk_size = len(power)/chunk_divisor
	for i in range(chunk_divisor):
		current_chunk = power[i*chunk_size:(i+1)*chunk_size] 
		power[i*chunk_size:(i+1)*chunk_size] = current_chunk-np.mean(current_chunk)
	return power

def plot_freq_data(power, xlabel = 'Frequency (Hz)', ylabel = 'Power (V$\cdot$s)'):
	FT = np.abs(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(power))))**2
	freqs = np.fft.fftshift(np.fft.fftfreq(len(power), 1))
	freq_res = round(1/len(power)*1000,2)
	print('Max power found at: ' + str(freqs[np.argmax(FT)]) + ' Hz')
	plt.figure()
	plt.plot(freqs, FT)
	plt.xlabel(xlabel, fontsize = 14)
	plt.ylabel(ylabel, fontsize = 14)
	plt.xlim(0,0.058)
	plt.title('Power Spectrum $\Delta$f = ' + str(freq_res) + ' kHz', fontsize = 18) 
	return freqs[np.argmax(FT)]


def get_max_freq(power, chunk_divisor):
	power = power[:int(len(power)/chunk_divisor)*chunk_divisor]
	chunk_size = len(power)/chunk_divisor
	max_freq = []
	times = []
	midpoint = chunk_divisor/2
	delta_freq = int(.001/(1./chunk_size))
	for i in range(chunk_divisor):
		times.append(ugradio.timing.lst(ugradio.timing.julian_date(time[i*chunk_size+midpoint]))-ra)
		if i == 0:
			FT = np.abs(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(power[i*chunk_size:(i+1)*chunk_size]))))**2
			freqs = np.fft.fftshift(np.fft.fftfreq(len(FT), 1))
			FT = FT[len(FT)/2:]
			freqs = freqs[len(freqs)/2:]
			max_freq.append(freqs[np.argmax(FT)])
		else:
			center_freq = int(max_freq[i-1]/(1./chunk_size))
			FT = np.abs(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(power[i*chunk_size:(i+1)*chunk_size]))))**2
			freqs = np.fft.fftshift(np.fft.fftfreq(len(FT), 1))
			FT = FT[len(FT)/2+center_freq-3*delta_freq:len(FT)/2+center_freq+3*delta_freq]
			freqs = freqs[len(freqs)/2+center_freq-3*delta_freq:len(freqs)/2+center_freq+3*delta_freq]
			max_freq.append(freqs[np.argmax(FT)])
	return times, max_freq

def guess_fringe(max_freqs, times):
	guess_freq = []
	guess = np.arange(.02, .045, 0.0001)
	midpoint = len(time)/len(max_freqs)/2
	S = [] 
	for j in range(len(guess)):
		ls = 0
		for i in range(len(max_freqs)):
			ls += (max_freqs[i]-guess[j]*np.cos(dec)*np.cos(times[i]))**2
		S.append(ls)
	return guess, S

def least_squares_2(measured, baseline, times):
	s_m = []
	t_m = []
	s_s = []
	t_t = []
	s_t = []
	s_y = []
	t_y = []

	
	for i in range(len(times)):
		times[i] = (ugradio.timing.julian_date(times[i]))
		sun_pos = ugradio.coord.sunpos(times[i])
		ra = sun_pos[0]*2*np.pi/360
		dec = sun_pos[1]*2*np.pi/360
		h_s = ugradio.timing.lst(times[i])-ra
		v_t = baseline[0]*np.cos(dec)*np.cos(times[i])+baseline[1]*np.cos(dec)*np.sin(ugradio.nch.lat)*np.sin(times[i])
		s_m.append(np.cos(2*np.pi*v_t))
		t_m.append(np.sin(2*np.pi*v_t))
		s_s.append(s_m[i]*s_m[i])
		s_t.append(s_m[i]*t_m[i])
		t_t.append(t_m[i]*t_m[i])
		s_y.append(s_m[i]*measured[i])
		t_y.append(t_m[i]*measured[i])

	s_m = np.sum(s_m)
	t_m = np.sum(t_m)
	s_s = np.sum(s_s)
	t_t = np.sum(t_t)
	s_t = np.sum(s_t)
	s_y = np.sum(s_y)
	t_y = np.sum(t_y)
	
	X = np.array([[s_s, s_t], [s_t, t_t]])
	X_T = X.transpose()
	Y = []
	Y.append(s_y)
	Y.append(t_y)
	Y = np.array(Y).transpose()
	
	alpha = np.dot(X_T, X)
	alpha_inv = np.linalg.inv(alpha)
	
	beta = np.dot(X_T, Y)

	a = np.dot(alpha_inv, beta)
	
	Y_BAR = np.dot(X, a)

	print("X")
	print(X)
	print("Y")
	print(Y)
	print("Y_BAR")
	print(Y_BAR)
	DEL_Y = Y-Y_BAR
	print(DEL_Y)

	s_2 = 1./(len(times)-2)*np.dot(DEL_Y.transpose(), DEL_Y)

	s_a_2 = []
	s_a_2.append(1*alpha_inv[0][0])
	s_a_2.append(1*alpha_inv[1][1])
	print(alpha_inv)
	print("a")
	print(a)
	print("s_a_2")
	print(s_a_2)
	return a

def least_squares(measured, times):
	s_m = []
	t_m = []
	s_s = []
	t_t = []
	s_t = []
	s_y = []
	t_y = []

	for i in range(len(times)):
		s_m.append(np.cos(dec)*np.cos(times[i]))
		t_m.append(np.cos(dec)*np.sin(ugradio.nch.lat)*np.sin(times[i]))
		s_s.append(s_m[i]*s_m[i])
		s_t.append(s_m[i]*t_m[i])
		t_t.append(t_m[i]*t_m[i])
		s_y.append(s_m[i]*measured[i])
		t_y.append(t_m[i]*measured[i])

	s_m = np.sum(s_m)
	t_m = np.sum(t_m)
	s_s = np.sum(s_s)
	t_t = np.sum(t_t)
	s_t = np.sum(s_t)
	s_y = np.sum(s_y)
	t_y = np.sum(t_y)
	
	X = np.array([[s_s, s_t], [s_t, t_t]])
	X_T = X.transpose()
	Y = []
	Y.append(s_y)
	Y.append(t_y)
	Y = np.array(Y).transpose()
	
	alpha = np.dot(X_T, X)
	alpha_inv = np.linalg.inv(alpha)
	
	beta = np.dot(X_T, Y)

	a = np.dot(alpha_inv, beta)
	
	Y_BAR = np.dot(X, a)

	print("X")
	print(X)
	print("Y")
	print(Y)
	print("Y_BAR")
	print(Y_BAR)
	DEL_Y = Y-Y_BAR
	print(DEL_Y)

	s_2 = 1./(len(times)-2)*np.dot(DEL_Y.transpose(), DEL_Y)

	s_a_2 = []
	s_a_2.append(1*alpha_inv[0][0])
	s_a_2.append(1*alpha_inv[1][1])
	print(alpha_inv)
	print("a")
	print(a)
	print("s_a_2")
	print(s_a_2)
	return a

def plot_guess(times, guess):
	guess_freqs = []
	for i in range(len(times)):
		guess_freqs.append(guess[0]*np.cos(dec)*np.cos(times[i])+guess[1]*np.cos(dec)*np.sin(times[i]*np.sin(ugradio.nch.lat)))
	plt.plot(times, guess_freqs)

data = remove_DC(data, chunk_divisor)
plt.plot(time, data)
plt.figure()
times, max_freq = get_max_freq(data, chunk_divisor/6)

max_freq = max_freq[:26]
times = times[:26]

plt.plot(times, max_freq, '--')
baseline = least_squares(max_freq, times)
print(baseline)
#baseline = [ 0.03822613, -0.02025612]
#plt.plot(time, data)
#plt.show()
#amplitudes = least_squares_2(data, baseline, time)

plot_guess(times, baseline)
#plt.plot(times, max_freq)

chunk_divisor = 165
time = time[0:chunk_divisor**2]
data = data[0:chunk_divisor**2]


def hamming_window(data, num_samples):
	for n in range(len(data)):
		current_window_value = 0.54-0.46*np.cos(2*np.pi*n/(num_samples-1))
		data[n] = current_window_value*data[n]
	return data

for i in range(len(data/chunk_divisor)):
	data[i*chunk_divisor:(i+1)*chunk_divisor] = hamming_window(data[i*chunk_divisor:(i+1)*chunk_divisor], chunk_divisor)

data = data.reshape((chunk_divisor, chunk_divisor))
plt.figure()
plt.imshow((np.fft.fftshift(np.abs(np.fft.fft(data))**2,axes=(1,))), cmap = 'gray_r')

sun = True
if sun == True:
	d = np.load('SunDataR51.npz')
	ra = 0
	chunk_divisor = 25

time = d['arr_0']
data = d['arr_1']

time = time[chunk_divisor**2:2*chunk_divisor**2]
data = data[chunk_divisor**2:2*chunk_divisor**2]

data = remove_DC(data, chunk_divisor)

#plt.plot(time, data)
plt.show()
amplitudes = least_squares_2(data, baseline, time)
def get_f_ideal(times, guesses):
	f_ideal = []
	for i in range(len(times)):
		sun_pos = ugradio.coord.sunpos(ugradio.timing.julian_date(times[i]))
		sun_ra = sun_pos[0]*2*np.pi/360
		sun_dec = sun_pos[1]*2*np.pi/360
		h_s = ugradio.timing.lst(ugradio.timing.julian_date(times[i]))-sun_ra
		f_ideal.append(guesses[0]*np.cos(sun_dec)*np.cos(h_s)+guesses[1]*np.cos(sun_dec)*np.sin(ugradio.nch.lat)*np.sin(h_s))

	return f_ideal
#times, max_freq = get_max_freq(data, chunk_divisor/5)
#plt.plot(times, max_freq)
f_ideal = get_f_ideal(time, baseline)
#plt.figure()
#plt.plot(f_ideal)
'''
plt.figure()
plt.plot(guess, least_squares)
plt.figure()
plt.plot(times[:26], max_freq)
plot_guess(times, guess[np.argmin(least_squares)])
plt.figure()
plt.plot(data)

data = data[:146*146]
data = data.reshape((146,146))
#data = data.reshape((chunk_divisor, chunk_divisor))
plt.figure()
plt.imshow(np.fft.fftshift(np.abs(np.fft.fft(data))**2,axes=(1,)))

plt.show()
'''

