from kinetics.brubach_2022 import BrubachModel
from utils.parameters import T0, P0, F_INLET
m = BrubachModel()
cov = m.solve_surface(T0, P0, F_INLET.copy())
print(cov)
