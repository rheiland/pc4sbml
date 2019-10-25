# Plot results in output file from:
# ./simulateSBML -t 10 -a -o toy.csv Toy_Model_for_PhysiCell.xml
#
import matplotlib.pyplot as plt
import numpy as np

#time,Energy,Glucose,Hydrogen,Oxygen,Alpha,Beta,Gamma,phenotype_cycle_transition_rate_0_1,phenotype_death_rate_0,E_cyc_threshold,E_death_threshold,ModelValue_5,ModelValue_6,Cell
#t,energy,glu,h,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt('toy.csv', delimiter=',', skiprows=1, unpack=True)
t,energy,glu,hydrogen,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt('species.csv', delimiter=',', skiprows=1, unpack=True)
plt.plot(t,energy, label='energy',color='black')
plt.plot(t,glu, label='glucose')
plt.plot(t,hydrogen, label='hydrogen')
plt.plot(t,o2, label='oxygen')

plt.xlabel('t')
#plt.ylabel('energy')
plt.title('Toy model')
plt.legend()
plt.show()
