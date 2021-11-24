import numpy as np
from sys import argv
from math import sqrt, sin, pi


def init_psi(n, N):
    xk = np.linspace(0, 1, N+1)
    psiI = np.zeros(N+1, dtype = np.float64)
    psiR = np.zeros(N+1, dtype = np.float64)
    for i, x in enumerate(xk):
        psiR[i] = sqrt(2) * sin(n*pi*x)

    print(psiR)
    print(psiI)
    return psiI, psiR


def H_operator(psi, N, k, w, tau):
    x = np.linspace(0, 1, N+1)
    H = np.zeros(N+1)
    for i in range(1, N):
        H[i] = -1./2 * (psi[i+1] + psi[i-1] -2*psi[i])/(1./N)**2 + k*(x[i] - 1./2)*psi[i]*sin(w * tau)

    return H

def simulate(psiR, psiI, N, k, w, S, dt):
    tau = 0
    x = np.linspace(0, 1, N+1)
    dx = 1./N
    with open(argv[2], 'w+') as out:
        out.write("t \t N \t x \t E \n")
        for s in range(S):
            psiR = psiR + H_operator(psiI, N, k, w, tau)*dt/2
            tau = tau + dt/2
            psiI = psiI - H_operator(psiR, N, k, w, tau)*dt
            tau = tau + dt/2
            psiR = psiR + H_operator(psiI, N, k, w, tau)*dt/2

            out.write("{:.3f}".format(tau) + "\t")
            out.write("{:.3f}".format(dx * np.sum(psiR**2 + psiI**2)) + "\t")
            out.write("{:.3f}".format(dx * np.sum(x * (np.square(psiR + np.square(psiI))))) + "\t")
            out.write("{:.3f}".format(dx * np.sum(psiR*H_operator(psiR, N, k, w, tau) + psiI*H_operator(psiI, N, k, w, tau))) + "\n")






if __name__ == "__main__":

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
    # HR = H_operator(psiR, int(params['N']), params['k'], params['w'], 0.0005 )
    # simulate(psiR, psiI, int(params['N']), params['k'], params['w'], int(params['S']), params['dt'])
