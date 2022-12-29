# Import tk module
import tkinter as tk
# Import themed tk
from tkinter import ttk
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)
from matplotlib.figure import Figure

icmstr="cm\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}"
degree="\N{DEGREE SIGN}"
kappastr="\u03BA"
Data=np.loadtxt('TD_Absorption.dat')
#Data[0,1]=Data[0,1]/2
icm2ifs=2.99792458e-5

def save():
#    np.savetxt('Absorption.txt',Ia)    
    update_canvas()

def load():
    Data=np.loadtxt('TD_Absorption.dat')

def quit():
    exit()

def update_kappa():
    kappa=0
#    kappa=1/(float(time.get())*float(sigma.get()))/3e-5
#    kappaLabel=tk.Label(text="Linebroadening parameter "+kappastr)
#    kappaValue=tk.Label(text=str(kappa))
#    kappaLabel.grid(row=12,column=0)
#    kappaValue.grid(row=12,column=1)

def lineshape():
    tc=float(time.get())
    si=float(sigma.get())
    t=np.linspace(0,tc*3,100)
    g=(si/tc)**2*(np.exp(-t/tc)+t/tc-1)
    I=np.fft.fft(g)

def update_canvas():
    # Plot Response function
    fig=Figure(figsize=(5,3),dpi=100)
    plot1 = fig.add_subplot(111)
    plot1.plot(Data[:,0],Data[:,1],'b')
    appf=np.exp(-Data[:,0]*float(homo.get())*icm2ifs)
    appf=appf*np.exp(-Data[:,0]**2*(float(inhomo.get())*icm2ifs)**2)
    plot1.plot(Data[:,0],Data[:,1]*appf,'r')
    # Draw arrows for dipoles
    plot1.set_ylabel("Response")
    plot1.set_xlabel("Time in fs")
# containing the Matplotlib figure
    fig.set_tight_layout(True)
    canvas = FigureCanvasTkAgg(fig,master = window)
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().grid(column=0,row=7,columnspan=2,sticky=tk.W+tk.E)

    # Do FFT
    Cdata=Data[:,1]+1j*Data[:,2]
    Cdata[0]=Cdata[0]/2
    I=np.fft.fft(Cdata)
    I=np.fft.fftshift(I)
    Ia=np.fft.fft(Cdata*appf)
    Ia=np.fft.fftshift(Ia)
    freq=np.fft.fftfreq(Data[:,0].shape[-1],Data[1,0])
    freq=np.fft.fftshift(freq)/icm2ifs

    cfreq=(float(minf.get())+float(maxf.get()))/2
    freq=freq+cfreq

    # Plot Spectrum
    fig=Figure(figsize=(5,3),dpi=100)
    plot1 = fig.add_subplot(111)
    plot1.plot(freq,I.real,'b')
    plot1.plot(freq,Ia.real,'r')
    # Draw arrows for dipoles
    plot1.set_ylabel("Intensity")
    plot1.set_xlabel("Wavenumbers in "+icmstr)
    plot1.set_xlim([float(minf.get()),float(maxf.get())])

# containing the Matplotlib figure
    fig.set_tight_layout(True)
    canvas = FigureCanvasTkAgg(fig,master = window)
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().grid(column=0,row=8,columnspan=2,sticky=tk.W+tk.E)
    # creating the Matplotlib toolbar
#    toolbar = NavigationToolbar2Tk(canvas,window)
#    toolbar.update()
    # placing the toolbar on the Tkinter window
#    canvas.get_tk_widget().grid()
    spec=open("Absorption.txt","w")
    for i in range(I.size):
        if (freq[i]>=float(minf.get()) and freq[i]<=float(maxf.get())):
            spec.write(str(freq[i])+" "+str(Ia[i].real)+" "+str(Ia[i].imag)+"\n")
    spec.close()

# Create and configure the window
window = tk.Tk()
# Add frame for options
#frame=tk.Frame(window)
# We want two columns with multiple rows
window.columnconfigure([0,1], minsize=150)
window.rowconfigure([0, 1,2,3,4,5,6,7,8,9,10,11,12,13], minsize=25)
window.title("Absorption Apodization")
window.resizable(width=False, height=False)

# Create Button to start Hamiltonian creation when pressed
btn_save = tk.Button(text="Save", command=save)
btn_update = tk.Button(text="Update", command=update_canvas) 
btn_load = tk.Button(text="Load", command=load)
btn_quit = tk.Button(text="Quit", command=quit)

# OptionMenu

# Average Frequency
homo=tk.Entry(window)
homo.insert(0, '0')
homoLabel=tk.Label(text="Homogeneous broadening in "+icmstr)
# Frequency Difference
inhomo=tk.Entry(window)
inhomo.insert(0, '0')
inhomoLabel=tk.Label(text="Inhomogeneous broadening in "+icmstr)
# Coupling
minf=tk.Entry(window)
minf.insert(0, '1100')
minfLabel=tk.Label(text="Min frequency in "+icmstr)
# Disorder magnitude
maxf=tk.Entry(window)
maxf.insert(0, '1300')
maxfLabel=tk.Label(text="Max frequency in "+icmstr)
# Correlation time
time=tk.Entry(window)
time.insert(0, '250')
timeLabel=tk.Label(text="Correlation time in fs")
# Correlation angle
angle=tk.Entry(window)
angle.insert(0, '90')
angleLabel=tk.Label(text="Correlation angle in "+degree)
# Timestep
deltat=tk.Entry(window)
deltat.insert(0, '10')
deltatLabel=tk.Label(text="Timestep in fs")
# Length
length=tk.Entry(window)
length.insert(0, '100000')
lengthLabel=tk.Label(text="Trajectory length in steps")
# Dipole angle
angle2=tk.Entry(window)
angle2.insert(0, '90')
angle2Label=tk.Label(text="Angle between dipoles in "+degree)
# Lineabroadening parameter
update_kappa()

# Configure Grid
homoLabel.grid(row=0,column=0)
homo.grid(row=0,column=1)
inhomoLabel.grid(row=1,column=0)
inhomo.grid(row=1,column=1)
minfLabel.grid(row=2,column=0)
minf.grid(row=2,column=1)
maxfLabel.grid(row=3,column=0)
maxf.grid(row=3,column=1)
#timeLabel.grid(row=4,column=0)
#time.grid(row=4,column=1)
#angleLabel.grid(row=5,column=0)
#angle.grid(row=5,column=1)
#deltatLabel.grid(row=6,column=0)
#deltat.grid(row=6,column=1)
#lengthLabel.grid(row=7,column=0)
#length.grid(row=7,column=1)
#angle2Label.grid(row=8,column=0)
#angle2.grid(row=8,column=1)
btn_update.grid(row=4,column=0)
#btn_load.grid(row=9,column=1)
btn_save.grid(row=5,column=0)
btn_quit.grid(row=5,column=1)
ttk.Separator(window,orient='horizontal').grid(row=6,column=0,columnspan=2,sticky='ew')
update_canvas()
window.mainloop()

