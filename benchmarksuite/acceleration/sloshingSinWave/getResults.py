#%%
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
import numpy as np
import seaborn as sns


caseStructure = [["gravity", "staticGravity","staticGravityPLIC","gravityRecon","gravityRecon2"],
                 ["Grid1", "Grid2", "Grid3", "Grid4", "Grid5"]]

baseCase = 'Cases'


solutionDir = 'maxU/0'
file = 'fieldMinMax.dat'

#%%
sol = casefoam.time_series(solutionDir,file,caseStructure,baseCase)

sol.reset_index(inplace=True)
sol.columns = ['t','minU', 'maxU','GravModel','Resolution']
sol = sol.drop(columns='minU')


maxU = sol.groupby(["GravModel","Resolution"]).max().reset_index()
maxU.to_csv("maxU.csv",index=False)
#%%
sns.lineplot(x="t",y="maxU",hue="GravModel",style="Resolution",data=sol)
plt.show()
# %%
