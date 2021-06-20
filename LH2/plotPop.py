import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('Pop.dat')

plt.plot(Data[:,0],Data[:,1])
plt.xlabel('Time [fs]',fontsize=16)
plt.ylabel('Population',fontsize=16)
plt.savefig('Population.png')
plt.show()

Data = np.loadtxt('PopF.dat')

s850=np.array([0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23,24,26])
s800=np.array([1, 4, 7, 10, 13, 16, 19, 22,25])
N=np.size(Data[:,0])
print(N)
P850=np.zeros(N)
for i in s850:
  index=1*27+i+1
  print(index)
  P850=P850+Data[:,index]

P800=np.zeros(N)
for i in s800:
  P800=P800+Data[:,1*27+i+1]

plt.plot(Data[:,0],P850)
plt.plot(Data[:,0],P800)
plt.xlabel('$\omega$ [cm$^{-1}$]',fontsize=16)
plt.ylabel('Population',fontsize=16)
plt.savefig('PopulationF.png')
plt.show()

