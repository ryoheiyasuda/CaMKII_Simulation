# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 08:10:29 2019

@author: yasudar
"""

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure #To use Figure command
from CaMKII_SimRev import CaMKII_Sim

(t1, activeCaMKII_WT, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(1, 1, 0.6, 0.1, False, False, 0)
(t1, activeCaMKII_T286A, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(0, 1, 0.6, 0.1, False, False, 0)
(t1, activeCaMKII_T286D, CaM_Bound, PhosphoCaMKII) = CaMKII_Sim(1, 0, 0.6, 0.1, False, False, 0)

#####################################################
#Plotting
#####################################################    

windowResolution = 100
plotsize = (8,4)
plotWindow = tk.Tk()
plotWindow.wm_title('Plot window')                
f = Figure(figsize = plotsize, dpi=windowResolution) #define the size of the figure.

a = f.add_subplot(1,2,1)
a.plot(t1, activeCaMKII_WT, '-k') #Black
a.plot(t1, activeCaMKII_T286D, '-m') #Magenta
a.plot(t1, activeCaMKII_T286A, '-g') #Green

a.set_xlim(-2, 20)
a.set_ylim(0, 0.8)

a = f.add_subplot(1,2,2)
a.plot(t1, activeCaMKII_WT, '-k') #Black
a.plot(t1, activeCaMKII_T286D, '-m') #Magenta
a.plot(t1, activeCaMKII_T286A, '-g') #Green

a.set_xlim(-10, 150)
a.set_ylim(0, 0.8)

canvas = FigureCanvasTkAgg(f, plotWindow)
canvas.get_tk_widget().pack(side = tk.TOP, fill = tk.BOTH, expand = True)

NavigationToolbar2Tk(canvas, plotWindow).update()
canvas._tkcanvas.pack(side = tk.TOP, fill = tk.BOTH, expand = True)
plotWindow.protocol("WM_DELETE_WINDOW", plotWindow.quit)
plotWindow.mainloop()
plotWindow.destroy()