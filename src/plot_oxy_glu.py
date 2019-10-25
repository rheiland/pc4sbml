# Plot results in output file from:
# ./simulateSBML -t 10 -a -o toy.csv Toy_Model_for_PhysiCell.xml
#
import sys
import matplotlib.pyplot as plt
import numpy as np

#fname = string.atof(sys.argv[1])
oname = 'oxygen.dat'
gname = 'glucose.dat'

#time,Energy,Glucose,Hydrogen,Oxygen,Alpha,Beta,Gamma,phenotype_cycle_transition_rate_0_1,phenotype_death_rate_0,E_cyc_threshold,E_death_threshold,ModelValue_5,ModelValue_6,Cell
#t,energy,glu,h,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt('toy.csv', delimiter=',', skiprows=1, unpack=True)

#t,energy= np.loadtxt('toy3.csv', delimiter=',', skiprows=1, unpack=True)
#t,energy= np.loadtxt('species.csv', delimiter=',', skiprows=0, unpack=True)
#t,spec= np.loadtxt('species.csv', delimiter=',', skiprows=0, unpack=True)
t,oxy= np.loadtxt(oname, delimiter=',', skiprows=0, unpack=True)
t,glu= np.loadtxt(gname, delimiter=',', skiprows=0, unpack=True)
#plt.plot(t,energy, label='energy',color='black')
plt.plot(t, oxy, label='oxygen', color='red')
plt.plot(t, glu, label='glucose', color='darkslateblue')  # https://matplotlib.org/examples/color/named_colors.html

plt.xlabel('t (sbmlsim)')
plt.title('sbmlsim results: (PhysiCell) t=0')
plt.legend()
plt.show()
