import matplotlib.pyplot as plt
import numpy as np

#file = input('Arquivo das médias: ')
file = './implicit1/meanFile.txt'
T, P = np.loadtxt(file, unpack=True, usecols=(0,1))
#T1, P1 = np.loadtxt('./implicit2/meanFile.txt', unpack=True, comments="#", usecols=(0,1))
#T2, P2 = np.loadtxt('./isoline/rho1.10/meanFile.txt', unpack=True)

# Transformar em Mega Pascal    
'''P = P*101325/(10**6)
#P1 = P1*101325/(10**6)
#P2 = P2*101325/(10**6)
'''
# Transformar em Kilo Bar
P = P*1.01325/1000
#P1 = P1*1.01325/1000


#plt.scatter(T1[:], P1[:], marker='^', color='blue', label=r'T 2000000 - $\rho$ = 1.05 g/cm^3')
#plt.plot(T1[:], P1[:], color='blue', linewidth=0.5)
#plt.scatter(T2, P2, marker='^', color='black', label=r'- $\rho$ = 1.10 g/cm^3')
#plt.plot(T2, P2, color='black', linewidth=0.5)
plt.scatter(T[:], P[:], marker='.', color='black', label=r'$\rho$ = 1.05 g/cm^3')
plt.plot(T[:], P[:], color='black', linewidth=0.5)
plt.ylabel('P (Kbar)')
plt.xlabel('T (K)')
plt.legend()
plt.title('TIP4Pe')
plt.ylim(bottom=0)
plt.show()