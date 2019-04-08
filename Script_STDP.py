# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 08:10:29 2019

@author: yasudar
"""

import tkinter as tk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure #To use Figure command
from CaMKII_SimRev import CaMKII_Sim

#Release probability = 0 or 1. Deterministic
(t1, activeCaMKII10, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(1, 1, 0.6, 0.1, False, True, 1)
(t1, activeCaMKII0, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(1, 1, 0.6, 0.1, False, True, 0)

#Under other situation, it is stochastic. We repeat it.
nRepeat = 20
siz = np.size(activeCaMKII10)
activeCaMKII5_all = np.zeros((siz, nRepeat))
activeCaMKII2_all = np.zeros((siz, nRepeat))

for i in range(nRepeat):
    (t1, activeCaMKII5, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(1, 1, 0.6, 0.1, False, True, 0.5)
    (t1, activeCaMKII2, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(1, 1, 0.6, 0.1, False, True, 0.2)
    activeCaMKII5_all[:,i] = activeCaMKII5
    activeCaMKII2_all[:,i] = activeCaMKII2

activeCaMKII5Ave = np.mean(activeCaMKII5_all, axis=1)
activeCaMKII2Ave = np.mean(activeCaMKII2_all, axis=1)


#####################################################
#Plotting
#####################################################    

windowResolution = 100
plotsize = (8,4)
plotWindow = tk.Tk()
plotWindow.wm_title('Plot window')                
f = Figure(figsize = plotsize, dpi=windowResolution) #define the size of the figure.

a = f.add_subplot(1,2,1)
a.plot(t1, activeCaMKII10, '-k')
a.plot(t1, activeCaMKII5_all, '-m')
a.plot(t1, activeCaMKII2_all, '-c')
a.plot(t1, activeCaMKII5Ave, '-r')
a.plot(t1, activeCaMKII2Ave, '-b')
a.plot(t1, activeCaMKII0, '-g')

a.set_xlim(-2, 20)
a.set_ylim(0, 0.8)

a = f.add_subplot(1,2,2)
a.plot(t1, activeCaMKII10, '-k')
a.plot(t1, activeCaMKII5_all, '-m')
a.plot(t1, activeCaMKII2_all, '-c')
a.plot(t1, activeCaMKII5Ave, '-r')
a.plot(t1, activeCaMKII2Ave, '-b')
a.plot(t1, activeCaMKII0, '-g')

a.set_xlim(-10, 150)
a.set_ylim(0, 0.8)

canvas = FigureCanvasTkAgg(f, plotWindow)
canvas.get_tk_widget().pack(side = tk.TOP, fill = tk.BOTH, expand = True)

NavigationToolbar2Tk(canvas, plotWindow).update()
canvas._tkcanvas.pack(side = tk.TOP, fill = tk.BOTH, expand = True)
plotWindow.protocol("WM_DELETE_WINDOW", plotWindow.quit)
plotWindow.mainloop()
plotWindow.destroy()