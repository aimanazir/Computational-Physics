# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 18:23:57 2020

@author: Ahmad Aiman Mohd Nazir 
"""
#import the relevent libraries
import matplotlib.pyplot as plt
from scipy.io.wavfile import read
from scipy.fftpack import rfft,irfft,rfftfreq
import numpy as np

#read the wav file of my favourite song
#the wav file in the form of mono
rate, data = read('TWICE.wav')

#calculate the time axis of the song
T = len(data)/rate
t = np.arange(0,T,1/rate)

#perform the real fast fourier transform and inverse real fast fourier transform
FT = rfft(data)
iFT = irfft(FT)

#calculate the frequency axis of song 
freq = rfftfreq(len(data),1/rate)

#plot the original waveform
plt.figure(1)
plt.plot(t,data)
plt.title('Original Waveform')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
#plot fast fourier transform 
plt.figure(2)
plt.plot(freq,FT)
plt.title('Fourier Transform')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
#plot inverse fourier transform
plt.figure(3)
plt.plot(t,iFT)
plt.title('Inverse Fourier Transform')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')

plt.show()
  

  