import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size']=8
plt.rcParams['ytick.major.size']=8
plt.rcParams['xtick.major.width']=2
plt.rcParams['ytick.major.width']=2
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'

Data = np.loadtxt('2D.par.dat')

w1Dim=int(np.sqrt(len(Data[:,0])))
w3Dim=w1Dim

wmin=11000
wmax=13000

step=500
mmax=np.max(np.abs(Data[:,3]))
Data[:,3]=Data[:,3]/mmax
Data[1,3]=1
Data[2,3]=-1

w1 = np.reshape(Data[:,0],(w1Dim,w3Dim))
w3 = np.reshape(Data[:,1],(w1Dim,w3Dim))

PlotMap = np.reshape(Data[:,3],(w1Dim,w3Dim))

fig=plt.figure()
plt.contourf(w1,w3,PlotMap,20,cmap=plt.cm.bwr)
plt.xlabel(r'$\omega_1$ [cm$^{-1}$]',fontsize=16)
plt.ylabel(r'$\omega_3$ [cm$^{-1}$]',fontsize=16)
plt.gca().set_aspect('equal')

plt.xlim(wmin,wmax)
plt.ylim(wmin,wmax)

plt.plot([ wmin , wmax],[wmin,wmax],color='black')
plt.xticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.yticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.show()
fig.savefig('2Dpar.png',dpi=400)

Data = np.loadtxt('2D.per.dat')

w1Dim=int(np.sqrt(len(Data[:,0])))
w3Dim=w1Dim

wmin=11000
wmax=13000

step=500
mmax=np.max(np.abs(Data[:,3]))
Data[:,3]=Data[:,3]/mmax
Data[1,3]=1
Data[2,3]=-1

w1 = np.reshape(Data[:,0],(w1Dim,w3Dim))
w3 = np.reshape(Data[:,1],(w1Dim,w3Dim))

PlotMap = np.reshape(Data[:,3],(w1Dim,w3Dim))

fig=plt.figure()
plt.contourf(w1,w3,PlotMap,20,cmap=plt.cm.bwr)
plt.xlabel(r'$\omega_1$ [cm$^{-1}$]',fontsize=16)
plt.ylabel(r'$\omega_3$ [cm$^{-1}$]',fontsize=16)
plt.gca().set_aspect('equal')

plt.xlim(wmin,wmax)
plt.ylim(wmin,wmax)

plt.plot([ wmin , wmax],[wmin,wmax],color='black')
plt.xticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.yticks(np.arange(wmin,wmax+1,step),fontsize=12,rotation=0)
plt.show()
fig.savefig('2Dper.png',dpi=400)
