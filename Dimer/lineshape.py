import numpy as np
import matplotlib.pyplot as plt

icm2ifs=3e-5
ifs2icm=1/icm2ifs
tc=500 # fs
si=50 # cm-1
w0=-500 # cm-1
N=1000
kappa=1/tc/si/icm2ifs
print(kappa)
if kappa<1:
    print('Slow modulation limit')
    tm=1/si/icm2ifs*10
if kappa>=1:    
    print('Fast modulation limit')
    tm=1/(si*tc*icm2ifs)**2*tc*10
t=np.linspace(0,tm,N)
dt=t[1]
gl=(si*tc*icm2ifs)**2*(np.exp(-t/tc)+t/tc-1)
g=np.exp(-1j*w0*icm2ifs*t)*np.exp(-gl)
I=np.fft.fft(g)[:N//2]
w=ifs2icm*np.fft.fftfreq(N,dt)[:N//2]*2*np.pi
plt.plot(w,np.real(I))
plt.plot(w,np.imag(I))
plt.show()
plt.plot(t,g)
plt.show()
