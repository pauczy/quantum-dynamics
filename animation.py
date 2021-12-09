#animacja dla przejścia rezonansowego w=9pi^2/2 dla n=4

import matplotlib

matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
from matplotlib import animation


sns.set_theme( style = 'whitegrid')

rho_rez = np.loadtxt('rho_4_9pi_n4.dat')
ffig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(-0.5, 3.6))
line, = ax.plot([], [], lw=2, c = 'midnightblue')
plt.xlabel(r'$\mathcal{x}$', fontsize=12, loc='right')
plt.ylabel(r'$\rho$', fontsize=12, loc='top')
plt.title(r'n = 4, $\kappa = 1$, $\omega = 9\pi^2/2$')

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

Writer = animation.writers['ffmpeg']
writer = Writer(fps=8, metadata=dict(artist='Me'), bitrate=-1, extra_args=['-vcodec', 'libx264'])
anim.save('resonance_n4.mp4', writer=writer)


plt.show()
