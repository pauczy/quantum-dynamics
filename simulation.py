import numpy as np
from math import sqrt, sin, pi
from numba import njit
from sys import argv
from time import time


def init_psi(n, N):
    xk = np.linspace(0, 1, N+1)
    psiI = np.zeros(N+1)
    psiR = np.zeros(N+1)
    for i in range(N+1):
        psiR[i] = sqrt(2) * sin(n*pi*xk[i])
    return psiI, psiR

@njit
def H_operator(psi, N, k, w, tau):

    dx = 1./N
    H = np.zeros(N+1)
    for i in range(1, N):
        H[i] = -1./2 * (psi[i+1] + psi[i-1] -2*psi[i])/dx**2 + k*(i/N - 1./2)*psi[i]*sin(w * tau)
    return H

def simulate(psiR, psiI, N, k, w, S, S_out, S_rho, dt):
    tau = 0
    x = np.linspace(0, 1, N+1)
    dx = 1./N
    with open(argv[2], 'w+') as out, open(argv[3], 'w+') as out_rho:
        out.write("t \t N \t x \t E \n")
        for s in range(S):
            psiR = psiR + H_operator(psiI, N, k, w, tau)*dt/2
            tau = tau + dt/2
            psiI = psiI - H_operator(psiR, N, k, w, tau)*dt
            tau = tau + dt/2
            psiR = psiR + H_operator(psiI, N, k, w, tau)*dt/2

            if (s % S_out == 0):
                out.write("{:.3f}".format(tau) + "\t")
                out.write("{:.3f}".format(dx * np.sum(psiR**2 + psiI**2)) + "\t")
                out.write("{:.3f}".format(dx * np.sum(x * (psiR**2 + psiI**2))) + "\t")
                out.write("{:.3f}".format(dx * np.sum(psiR*H_operator(psiR, N, k, w, tau) + psiI*H_operator(psiI, N, k, w, tau))) + "\n")

            if (s % S_rho == 0):
                rho = psiR[0:N+1:2]**2 + psiI[0:N+1:2]**2
                np.savetxt(out_rho, rho, delimiter = '\t')



if __name__ == "__main__":

    tic = time()

    params = {}
    with open(argv[1]) as f:
        for line in f:
            val, key = line.split("#")
            params[key.strip()] = float(val.strip())

    psiI, psiR = init_psi(int(params['n']), int(params['N']))
    xk = np.linspace(0, 1, int(params['N'])+1)
    N = 1./int(params['N']) * np.sum(psiR**2 + psiI**2)
    x = 1./int(params['N']) * np.sum(xk * (psiR**2 + psiI**2))
    E = 1./int(params['N']) + np.sum(psiR * H_operator(psiR, int(params['N']), params['k'], params['w'], 0) + psiI * H_operator(psiI, int(params['N']), params['k'], params['w'], 0))
    print(N)
    print(x)
    print(E)
    simulate(psiR, psiI, int(params['N']), params['k'], params['w'], int(params['S']), int(params['S_out']), int(params['S_rho']), params['dt'])

    toc = time()
    print('execution time: ', (toc - tic))
