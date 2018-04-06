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

def calc_mf_theory(ff, time):
	#R_f = np.radians(np.linspace(.45, .55, len(time)))
	R_f = np.linspace(0,2.5, len(time))
	ff = np.array(ff)
	N = 200
	vals = np.linspace(-N,N,2*N+1)
	mf = 0
	total = 0
#	R = np.radians(.5)
	
	for n in range(len(vals)):
		mf+=10**(-4)*np.sqrt(1-(vals[n]/N)**2)*np.cos(2*np.pi*R_f*vals[n]/N)	
	plt.plot(R_f, mf)
	plt.figure()
	return mf

def get_ff(times, coeff):
	ff = []
	dec = 0
	for i in range(len(times)):
		ff.append(coeff[0]*(np.cos(dec)*np.cos(times[i])) + coeff[1]*(np.cos(dec)*np.sin(ugradio.nch.lat)*np.sin(times[i])))
	return ff

time = np.load('sun_ha.npz')['arr_1']
coeff = [ 0.03822613, 0]
ff = get_ff(time, coeff)
mf = calc_mf_theory(ff, time)
data = remove_DC(data, chunk_divisor)
plt.plot(time, ff)
plt.figure()
plt.plot(ff, data[::-1])
plt.figure()
plt.plot(time, mf)
plt.show()
