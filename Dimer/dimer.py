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
icmstr="cm-1"
#degree="o"
#kappastr="k"

def save():
    energy=open("Energy.txt","w")
    dipole=open("Dipole.txt","w")
    summary=open("Summary.txt","w")
    # Find correlation parameters
    a=np.exp(-float(deltat.get())/float(time.get()))
    b=np.sqrt(1-a*a)
    sig=float(sigma.get())
    JJ=J.get()
    cang=np.cos(float(angle.get())*np.pi/180)
    sang=np.sqrt(1-cang*cang)
    cang2=np.cos(float(angle2.get())*np.pi/180)
    sang2=np.sqrt(1-cang2*cang2)
    wao=random.gauss(0,sig)
    wbo=random.gauss(0,sig)
    for i in range(int(length.get())):
        wa=wao*a+b*random.gauss(0,sig)
        wb=wbo*a+b*random.gauss(0,sig)
        w1=wa-float(deltaw.get())/2+float(w0.get())
        w2=wa*cang+wb*sang+float(deltaw.get())/2+float(w0.get())
        energy.write(str(i)+" "+str(w1)+" "+JJ+" "+str(w2)+"\n")
        dipole.write(str(i)+" 1.0 "+str(cang2)+" 0.0 "+str(sang2)+" 0.0 0.0\n")
        wao=wa
        wbo=wb      
    energy.close()
    dipole.close()
    summary.write("Average energy: "+str(w0.get())+"\n")
    summary.write("Frequency difference: "+str(deltaw.get())+"\n")
    summary.write("Coupling: "+str(J.get())+"\n")
    summary.write("Disorder magnitude: "+str(sigma.get())+"\n")
    summary.write("Correlation time: "+str(time.get())+"\n")
    summary.write("Correlation angle: "+str(angle.get())+"\n")
    summary.write("Timestep: "+str(deltat.get())+"\n")
    summary.write("Trajectory length: "+str(length.get())+"\n")
    summary.write("Angle between dipoles: "+str(angle2.get())+"\n")
    summary.close()
    update_canvas()

def load():
    summary=open("Summary.txt","r")
    line=summary.readline()
    w0.delete(0,tk.END)
    w0.insert(0,line.split(' ')[2])
    line=summary.readline()
    deltaw.delete(0,tk.END)
    deltaw.insert(0,line.split(' ')[2])
    line=summary.readline()
    J.delete(0,tk.END)
    J.insert(0,line.split(' ')[1])
    line=summary.readline()
    sigma.delete(0,tk.END)
    sigma.insert(0,line.split(' ')[2])
    line=summary.readline()
    time.delete(0,tk.END)
    time.insert(0,line.split(' ')[2])
    line=summary.readline()
    angle.delete(0,tk.END)
    angle.insert(0,line.split(' ')[2])
    line=summary.readline()
    deltat.delete(0,tk.END)
    deltat.insert(0,line.split(' ')[1])
    line=summary.readline()
    length.delete(0,tk.END)
    length.insert(0,line.split(' ')[2])
    line=summary.readline()
    angle2.delete(0,tk.END)
    angle2.insert(0,line.split(' ')[3])
    summary.close()

def quit():
    exit()

def update_kappa():
    kappa=1/(float(time.get())*float(sigma.get()))/3e-5
    kappaLabel=tk.Label(text="Linebroadening parameter "+kappastr)
    kappaValue=tk.Label(text=str(kappa))
    kappaLabel.grid(row=12,column=0)
    kappaValue.grid(row=12,column=1)

def lineshape():
    tc=float(time.get())
    si=float(sigma.get())
    t=np.linspace(0,tc*3,100)
    g=(si/tc)**2*(np.exp(-t/tc)+t/tc-1)
    I=np.fft.fft(g)

def update_canvas():
    update_kappa()
    fig=Figure(figsize=(5,5),dpi=100)
    plot1 = fig.add_subplot(111)
    xmin=float(w0.get())-float(J.get())-3*float(sigma.get())-0.5*float(deltaw.get())
    xmax=2*float(w0.get())-xmin
    x = np.linspace(xmin,xmax,500)
    y1 = np.exp(-(x-float(w0.get())-float(deltaw.get())/2)**2/2/float(sigma.get())**2)
    y2 = np.exp(-(x-float(w0.get())+float(deltaw.get())/2)**2/2/float(sigma.get())**2)
    plot1.plot(y1,x,'b')
    plot1.plot(y2,x,'r')
    # Draw arrow for coupling
    dw=float(deltaw.get())
    if dw<float(sigma.get())*2:
        dw=float(sigma.get())
    a1=float(w0.get())+dw/2
    a2=a1-dw
    plot1.annotate(text='',xy=(0.5,a1),xytext=(0.5,a2),arrowprops=dict(arrowstyle='<->'))
    Jtxt="J="+J.get()+" "+icmstr
    plot1.annotate(text=Jtxt,xy=(0.6,(a1+a2)/2),xytext=(0.6,(a1+a2)/2))
    # Draw arrows for dipoles
    ###

    plot1.set_ylabel("Wavenumber in "+icmstr)
# containing the Matplotlib figure
    fig.set_tight_layout(True)
    canvas = FigureCanvasTkAgg(fig,master = window)
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().grid(column=0,row=13,columnspan=2,sticky=tk.W+tk.E)

    # creating the Matplotlib toolbar
#    toolbar = NavigationToolbar2Tk(canvas,window)
#    toolbar.update()
    # placing the toolbar on the Tkinter window
#    canvas.get_tk_widget().grid()

# Create and configure the window
window = tk.Tk()
# Add frame for options
#frame=tk.Frame(window)
# We want two columns with multiple rows
window.columnconfigure([0,1], minsize=150)
window.rowconfigure([0, 1,2,3,4,5,6,7,8,9,10,11,12,13], minsize=25)
window.title("Dimer model creator")
window.resizable(width=False, height=False)

# Create Button to start Hamiltonian creation when pressed
btn_save = tk.Button(text="Save", command=save)
btn_update = tk.Button(text="Update", command=update_canvas) 
btn_load = tk.Button(text="Load", command=load)
btn_quit = tk.Button(text="Quit", command=quit)

# OptionMenu

# Average Frequency
w0=tk.Entry(window)
w0.insert(0, '1200')
w0Label=tk.Label(text="Average frequeny in "+icmstr)
# Frequency Difference
deltaw=tk.Entry(window)
deltaw.insert(0, '0')
deltawLabel=tk.Label(text="Frequeny difference in "+icmstr)
# Coupling
J=tk.Entry(window)
J.insert(0, '25')
JLabel=tk.Label(text="Coupling in "+icmstr)
# Disorder magnitude
sigma=tk.Entry(window)
sigma.insert(0, '25')
sigmaLabel=tk.Label(text="Disorder magnitude in "+icmstr)
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
w0Label.grid(row=0,column=0)
w0.grid(row=0,column=1)
deltawLabel.grid(row=1,column=0)
deltaw.grid(row=1,column=1)
JLabel.grid(row=2,column=0)
J.grid(row=2,column=1)
sigmaLabel.grid(row=3,column=0)
sigma.grid(row=3,column=1)
timeLabel.grid(row=4,column=0)
time.grid(row=4,column=1)
angleLabel.grid(row=5,column=0)
angle.grid(row=5,column=1)
deltatLabel.grid(row=6,column=0)
deltat.grid(row=6,column=1)
lengthLabel.grid(row=7,column=0)
length.grid(row=7,column=1)
angle2Label.grid(row=8,column=0)
angle2.grid(row=8,column=1)
btn_update.grid(row=9,column=0)
btn_load.grid(row=9,column=1)
btn_save.grid(row=10,column=0)
btn_quit.grid(row=10,column=1)
ttk.Separator(window,orient='horizontal').grid(row=11,column=0,columnspan=2,sticky='ew')
update_canvas()
window.mainloop()

