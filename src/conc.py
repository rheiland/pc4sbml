import sys
import matplotlib.pyplot as plt
import numpy as np

#fname = sys.argv[1]
#print('csv_file=',fname)

#time,Energy,Glucose,Hydrogen,Oxygen,Alpha,Beta,Gamma,phenotype_cycle_transition_rate_0_1,phenotype_death_rate_0,E_cyc_threshold,E_death_threshold,ModelValue_5,ModelValue_6,Cell
#t,energy,glu,h,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt('toy.csv', delimiter=',', skiprows=1, unpack=True)
t,oxygen = np.loadtxt("oxygen.dat", delimiter=',', skiprows=0, unpack=True)
t,glucose = np.loadtxt("glucose.dat", delimiter=',', skiprows=0, unpack=True)
t,energy = np.loadtxt("energy.dat", delimiter=',', skiprows=0, unpack=True)

plt.plot(t,oxygen, label='oxygen')
plt.plot(t,glucose, label='glucose')
plt.plot(t,energy, label='energy',color='black')

plt.xlabel('t')
#plt.ylabel('energy')
#plt.title('Toy model')
plt.title("from solve_sbml, hetero, cellID=0, t=0")
plt.legend()
#png_fname = fname[:-4]+'.png'
png_fname = 'oge.png'
print(png_fname)
plt.savefig(png_fname)
plt.show()
