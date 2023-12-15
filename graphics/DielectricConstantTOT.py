import numpy as np
import matplotlib.pyplot as plt

def plot_DielectricTemperature(file):
    
    T, e = np.loadtxt(file, unpack=True)
    
    eref = [90, 79, 73]
    
    plt.scatter(T, e, label=r'$\rho$ = 1.05 g/cm^3')
    plt.plot(T, e, linewidth=0.5)
    plt.scatter([280, 300, 320], eref, color='black', label=r'$\rho_{ref}$ = 1.00 g/cm^3')
    plt.plot([280, 300, 320], eref, color='black', linewidth=0.5)
    plt.ylabel(r'$\epsilon$')
    plt.xlabel('T (K)')
    plt.title('Contante Dielétrica')
    #plt.ylim(bottom=30, top=130)
    #plt.xlim(left=220, right=330)
    plt.legend()
    plt.show()
    
def read_DipoleData(filename):

    mux = []
    muy = []
    muz = []

    with open(filename, 'r') as f:

        lines = f.readlines()[4:]

        for i in range(len(lines)-3):

            mux.append(float(lines[i].split(' ')[1])*angs*1.6*(10**-19))
            muy.append(float(lines[i+1].split(' ')[1])*angs*1.6*(10**-19))
            muz.append(float(lines[i+2].split(' ')[1])*angs*1.6*(10**-19))

            i = i + 4

        return mux, muy, muz


def plot_Dipole(mux, muy, muz, T):

    mu = []
    time = []
    for i in range(len(mux)):
        mu.append((mux[i]**2 + muy[i]**2 + muz[i]**2)**(1/2))
        time.append(i)

    plt.plot(time, mu)
    plt.xlabel('Time (fs)')
    plt.ylabel('M (C.m)')
    plt.title(f'Momento de dipolo total (T = {T}K)')
    plt.show()

    return mu

def dielectric_constante(mu, V, T):

    kB = 1.38*(10**(-23))
    e0 = 8.85*(10**(-12))

    M2 = 0
    Mmean = 0
    for m in mu:
        M2 = M2 + m*m
        Mmean = Mmean + m
        
    time = len(mu)

    Mvar = (M2/time) - (Mmean/time)**2
    
    #e = 1 + (Mvar)/(3*e0*kB*V*T)
    #e = 1 + (4*np.pi*Mvar)/(3*kB*V*T)
    e = 1 + (4*np.pi*Mvar)/(3*e0*kB*V*T)
    
    print(f'{T} {e}')
    
def graph_M():

    mu = np.loadtxt("./implicit/dipole.txt", unpack=True)

    time = np.arange(len(mu))

    plt.plot(time, mu)
    plt.xlabel("Time")
    plt.ylabel("M")
    plt.title("Momento de Dipolo")
    plt.show()
    
    return mu

def save_M(mu):

    with open("./implicit/mu.txt", "a") as f:
        
        for m in mu:
            f.write(f"{m}\n")
        
    f.close()

global angs

angs = 10**(-10)
lbox = 14.2*angs # m
V = lbox**3
T = 320
i = 0

mux, muy, muz = read_DipoleData(f'../../TIP4Pe/implicit/2dipoleMomentT{T}.sh')

#mu = plot_Dipole(mux[i:], muy[i:], muz[i:], T)

#save_M(mu)

#dielectric_constante(mu, V, T)

#tot_mu = graph_M()

#dielectric_constante(tot_mu, V, T)


plot_DielectricTemperature('./implicit/dielectric.txt')



