# Plot results in output file from:
# ./simulateSBML -t 10 -a -o toy.csv Toy_Model_for_PhysiCell.xml
#
import sys
import matplotlib.pyplot as plt
import numpy as np

fname = sys.argv[1]
print('csv_file=',fname)

#time,Energy,Glucose,Hydrogen,Oxygen,Alpha,Beta,Gamma,phenotype_cycle_transition_rate_0_1,phenotype_death_rate_0,E_cyc_threshold,E_death_threshold,ModelValue_5,ModelValue_6,Cell
#t,energy,glu,h,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt('toy.csv', delimiter=',', skiprows=1, unpack=True)
t,energy,glu,hydrogen,o2,alpha,beta,gamma,rate01,death_rate0,E_cyc_thresh,E_death_thresh,val5,val6,cell = np.loadtxt(fname, delimiter=',', skiprows=1, unpack=True)
plt.plot(t,o2, label='oxygen')
plt.plot(t,glu, label='glucose')
# plt.plot(t,hydrogen, label='hydrogen')
plt.plot(t,energy, label='energy',color='black')

plt.xlabel('t')
#plt.ylabel('energy')
#plt.title('Toy model')
plt.title(fname)
plt.legend()
png_fname = fname[:-4]+'.png'
print(png_fname)
plt.savefig(png_fname)
plt.show()
