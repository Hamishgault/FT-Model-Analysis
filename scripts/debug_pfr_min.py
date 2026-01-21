from kinetics.brubach_2022 import BrubachModel
from reactor.pfr import PFR
from utils.parameters import T0,P0
from utils.species import SPECIES_IDX
import numpy as np
m=BrubachModel()
pfr=PFR(1e-3,1.0,m,nu=None)
F0=np.zeros(len(SPECIES_IDX))
F0[SPECIES_IDX['CO']]=1.0
F0[SPECIES_IDX['H2']]=3.0
F0[SPECIES_IDX['CO2']]=0.1
z_eval=np.linspace(0,1,11)
sol=pfr.run(T0,P0,F0,z_eval=z_eval)
print('min,max overall', sol.y.min(), sol.y.max())
print('\nper species mins:')
names=['CO','H2','CH4','C2_4','C5plus','H2O','CO2']
for i in range(len(names)):
    print(names[i], sol.y[i,:].min())
print('\nAny negatives? ', (sol.y < -1e-12).any())
print('min negative value', sol.y[sol.y < -1e-12].min() if (sol.y < -1e-12).any() else None)