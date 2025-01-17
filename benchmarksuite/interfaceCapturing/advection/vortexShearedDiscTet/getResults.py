# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns

sns.set_style("ticks")

caseStructure = [['isoAlpha', 'plicRDF','MULES'],
                 ['x32', 'x64', 'x128', 'x256', 'x512']]

baseCase = 'Cases'
solutionDir = 'volumeFractionError/0'
file = 'volumeFractionError.dat'

error = casefoam.time_series(solutionDir,file,caseStructure,baseCase)
error = error.reset_index("t")

error = error[error["t"] == 8]
error = error.rename(columns={"var_0": "Model", "var_1": "Res"})
error = error.rename(columns={1: "E1", 2: "Emass", 3: "Ebound", 4 : "delta"}) # delta = avg cell dist

error["nCells"] = 1/error["delta"] # domain size is 1
error['Res'] = error['Res'].str.replace('x', '')
error['Res'] = error['Res'].apply(pd.to_numeric)

print("error",error)
error.to_csv("error.csv",index=False)
sns.lineplot(x="Res", y="E1",markers=True,hue="Model",style="Model",data=error)
plt.xscale('log')
plt.yscale('log')
plt.savefig("error.png")
# plt.show()

# save profiling data

profData = casefoam.profiling(8,processorDir="processor0",caseStructure=caseStructure,baseCase=baseCase)
profData = profData.rename(columns={"var_0":"Method","var_1":"Res"})
mask = profData["description"].str.contains('::advect\(')
profData.loc[mask,["Method","Resolution","totalTime"]].to_csv("profData.csv",index=False)

