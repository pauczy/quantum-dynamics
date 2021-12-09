#animacja dla przejścia rezonansowego w=9pi^2/2 dla n=4

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from matplotlib import animation


sns.set_theme( style = 'whitegrid')

rho_rez = np.loadtxt('rho_4_9pi_n4.dat')
ffig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-0.5, 3.5))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(0, 1, 51)
    y = rho_rez[51*i:51*(i+1)]
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(ffig, animate, init_func=init,
                               frames=1000, interval=200, blit=True)
plt.show()
