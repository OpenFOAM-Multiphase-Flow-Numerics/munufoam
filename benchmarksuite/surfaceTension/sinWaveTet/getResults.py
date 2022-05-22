# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns
from sinwave_prosperetti import sinwave_prosperetti


cases = [['geoVoF_Brackbill','algVoF_Brackbill'],
        ['coarse', 'mid']]#, 'fine']]

baseCase = 'Cases'
solutionDir = 'surfaces'
file = 'alpha.water_constantIso.raw'

rho1 = 1  # Density liquid 1
rho2 = 1  # Density liquid 2
wavelength = 0.003  # Wavelength
H0 = 3e-5  # Initial height
nu = 0.001  # Kinematic viscosity
sigma = 1       # Surface tension
ana = sinwave_prosperetti(rho1, rho2, wavelength, H0, nu, sigma)
analytical = ana.tabulatedData(duration=25, N=101,)
analytical.columns = ['time', "analytical"]
analytical["analytical"] = abs(analytical["analytical"])


postFunction = postFunctions.getFreeSurfaceWallAndCentre

sol = casefoam.posField_to_timeSeries(
    solutionDir, file, postFunction, cases, baseCase, axis=1)
sol = sol.reset_index()
sol.columns = ['time', 'min', 'mean', 'max',
               'Method', 'nCells']
sol = sol.replace('coarse', 32)
sol = sol.replace('mid', 64)
sol = sol.replace('fine', 128)
sol['max'] /= 3e-5
sol['time'] /= 1.475e-5

ax = analytical.plot(style='.', x='time', y='analytical',
                     c='black', marker='+', ms=7)
sns.set_style("ticks")
ax = sns.lineplot(x='time', y='max', hue='Method', style="nCells", data=sol, ax=ax)
plt.ylabel('Relative amplitude')
plt.xlabel('Non-dimensional time')
plt.savefig("sinWaveTet.pdf")
sol.to_csv("sinWaveTet.csv", index=False)



plt.show()
