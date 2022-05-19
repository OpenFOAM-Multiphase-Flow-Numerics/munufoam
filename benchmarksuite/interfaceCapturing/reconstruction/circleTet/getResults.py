# -*- coding: utf-8 -*-
from black import err
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
import seaborn as sns

sns.set_style("ticks")

caseStructure = [
    ["isoAlpha", "plicRDF", "gradAlpha"],
    ['x8','x16','x32', 'x64', 'x128', 'x256', 'x512', 'x1024'],
]

baseCase = "Cases"
solutionDir = "reconstructionError/0"
file = "error.dat"

error = casefoam.time_series(solutionDir, file, caseStructure, baseCase)
error = error.reset_index(drop=True)

error = error.rename(
    columns={
        1: "LNormal_1",
        2: "LNormal_Inf",
        3: "LCentre_1",
        4: "LCentre_Inf",
        9: "delta", # delta = avg cell dist
        10: "nCells",
        "var_0": "Model",
        "var_1": "Res",
    }
)  
error = error.drop([5,6,7,8,11],axis=1)

error['Res'] = error['Res'].str.replace('x', '')
error['Res'] = error['Res'].apply(pd.to_numeric)

error.to_csv("error.csv",index=False)
sns.lineplot(x="Res", y="LNormal_Inf",markers=True,hue="Model",style="Model",data=error)
plt.xscale('log')
plt.yscale('log')
plt.savefig("LNormal_Inf.png")

plt.figure()
sns.lineplot(x="Res", y="LNormal_1",markers=True,hue="Model",style="Model",data=error)
plt.xscale('log')
plt.yscale('log')
plt.savefig("LNormal_1.png")

plt.figure()
sns.lineplot(x="Res", y="LCentre_Inf",markers=True,hue="Model",style="Model",data=error)
plt.xscale('log')
plt.yscale('log')
plt.savefig("LCentre_Inf.png")

plt.figure()
sns.lineplot(x="Res", y="LCentre_1",markers=True,hue="Model",style="Model",data=error)
plt.xscale('log')
plt.yscale('log')
plt.savefig("LCentre_1.png")

# # save profiling data

profData = casefoam.profiling(0.001,caseStructure=caseStructure,baseCase=baseCase)
profData = profData.rename(columns={"var_0":"Method","var_1":"Resolution"})
mask = profData["description"].str.contains('::reconstruct\(')
profData = profData[mask]
profData = profData[profData["calls"] > 1] # reconstruct is also called in the constructor
profData.loc[:,["Method","Resolution","totalTime"]].to_csv("profData.csv",index=False)
