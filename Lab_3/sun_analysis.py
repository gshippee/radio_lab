import ugradio
import matplotlib.pylab as plt
import numpy as np
import argparse
import os

sun = True
if sun == True:
	d = np.load('SunDataR51.npz')
	ra = 0
	chunk_divisor = 144

time = d['arr_0']
data = d['arr_1']
print(len(data))

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

def plot_guess(times, guess):
	guess_freqs = []
	for i in range(len(times)):
		guess_freqs.append(guess[0]*np.cos(dec)*np.cos(times[i])+guess[1]*np.cos(dec)*np.sin(times[i]*np.sin(ugradio.nch.lat)))
	plt.plot(times, guess_freqs)

def calc_mf_theory(ff):
	R = np.radians(np.linspace(0.45, 0.55, .1))
	N = 200
	n = np.linspace(-N,N,2*N+1)
	mf = 0
	total = 0
	for j in range(len(R)):
		for i in range(len(n)):
			mf+=np.sqrt(1-(n/N)**2)*np.cos(2*np.pi*ff*R[j]*n/N)	
		plt.plot(mf)


def get_ff(times, coeff):
	ff = []
	ff.append(cooeff[0]*(np.cos(dec)*np.cos(times[i])) + coeff[1]*(np.cos(dec)*np.sin(ugradio.nch.lat)*np.sin(times[i])))
	

def fit_squares(h_a, dec, max_freqs):
	s_m = []
	t_m = []
	s_s = []
	t_t = []
	s_t = []
	s_y = []
	t_y = []

	for i in range(len(h_a)):
		dec[i] = 0
		s_m.append(np.cos(dec[i])*np.cos(h_a[i]))
		t_m.append(np.cos(dec[i])*np.sin(ugradio.nch.lat)*np.sin(h_a[i]))
		s_s.append(s_m[i]*s_m[i])
		s_t.append(s_m[i]*t_m[i])
		t_t.append(t_m[i]*t_m[i])
		s_y.append(s_m[i]*max_freqs[i])
		t_y.append(t_m[i]*max_freqs[i])

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
	a = [14.5,0]
	
	Y_BAR = np.dot(X, a)

	print("X")
	print(X)
	print("Y")
	print(Y)
	print("Y_BAR")
	print(Y_BAR)
	DEL_Y = Y-Y_BAR
	print(DEL_Y)
	print(DEL_Y)

	s_2 = 1./(len(h_a)-2)*np.dot(DEL_Y.transpose(), DEL_Y)

	s_a_2 = []
	s_a_2.append(1*alpha_inv[0][0])
	s_a_2.append(1*alpha_inv[1][1])
	print(alpha_inv)
	print("a")
	print(a)
	print("s_a_2")
	print(s_a_2)
	return a


def get_sun_ha(time):
	h_a = []
	dec = []
	for i in range(len(time)):
		time[i] = (ugradio.timing.julian_date(time[i]))
		sun_pos = ugradio.coord.sunpos(time[i])
		ra = sun_pos[0]*2*np.pi/360
		dec.append(sun_pos[1]*2*np.pi/360)
		h_s = ugradio.timing.lst(time[i])-ra
		h_a.append(h_s)
	np.savez('sun_ha', h_s, dec)
	

def get_max_freq(power, time, chunk_divisor):
	power = power[:int(len(power)/chunk_divisor)*chunk_divisor]
	chunk_size = len(power)/chunk_divisor
	max_freq = []
	h_a = []
	dec = []
	midpoint = chunk_divisor/2
	print(chunk_size)
	delta_freq = int(.001/(1./chunk_size))
	print(delta_freq)
	for i in range(chunk_divisor):
		time[i] = (ugradio.timing.julian_date(time[i]))
		sun_pos = ugradio.coord.sunpos(time[i])
		ra = sun_pos[0]*2*np.pi/360
		dec.append(sun_pos[1]*2*np.pi/360)
		h_s = ugradio.timing.lst(time[i])-ra
		h_a.append(h_s)
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
		print(freqs[np.argmax(FT)])
	return h_a, dec, max_freq


def find_coeff(h_a, dec, data):
	s_m = []
	t_m = []
	s_s = []
	t_t = []
	s_t = []
	s_y = []
	t_y = []
	for i in range(len(h_a)):
		s_m.append(np.cos(dec[i])*np.cos(h_a[i]))
		t_m.append(np.cos(dec[i])*np.sin(ugradio.nch.lat)*np.sin(h_a[i]))
	Y = []
	Y.append(data)
	Y = np.array(Y).transpose()

	X = []
	X.append(s_m)
	X.append(t_m)
	X_T = np.array(X)
	X = X_T.transpose()
	alpha = np.linalg.inv(np.dot(X_T, X))
	print(alpha.shape)
	beta = np.dot(X_T,Y)
	print(beta.shape)	
	coeff = np.dot(alpha, beta) 
	print(coeff)

def least_squares(actual, theor):
	s = 0
	for i in range(len(actual)):
		s+=((actual[i]-theor[i])**2)
	return s

R = 1.3914*10**9

data = remove_DC(data, chunk_divisor)

#get_sun_ha(time)
h_a, dec, max_freq = get_max_freq(data, time, chunk_divisor/8)


#calc_mf_theory(max_freq, R, 2000)

guess=fit_squares(h_a, dec, max_freq)
find_coeff(h_a,dec,max_freq)
plt.plot(h_a, max_freq)
#plot_guess(h_a, guess)
#plt.figure()
#plt.plot(time, data)

#plot_freq_data(data)

plt.show()
