# Plot results in output file from:
# ./simulateSBML -t 10 -a -o toy.csv Toy_Model_for_PhysiCell.xml
#
import sys
import matplotlib.pyplot as plt
import numpy as np

#fname = string.atof(sys.argv[1])
fname = sys.argv[1]

#time,Energy,Glucose,Hydrogen,Oxygen,Alpha,Beta,Gamma,phenotype_cycle_transition_rate_0_1,phenotype_death_rate_0,E_cyc_threshold,E_death_threshold,ModelValue_5,ModelValue_6,Cell
#t,energy,glu,h,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt('toy.csv', delimiter=',', skiprows=1, unpack=True)

#t,energy= np.loadtxt('toy3.csv', delimiter=',', skiprows=1, unpack=True)
#t,energy= np.loadtxt('species.csv', delimiter=',', skiprows=0, unpack=True)
#t,spec= np.loadtxt('species.csv', delimiter=',', skiprows=0, unpack=True)
t,spec= np.loadtxt(fname, delimiter=',', skiprows=0, unpack=True)
#plt.plot(t,energy, label='energy',color='black')
plt.plot(t, spec, label=fname,color='black')

plt.xlabel('t')
plt.title('from sbmlsim')
plt.legend()
plt.show()
