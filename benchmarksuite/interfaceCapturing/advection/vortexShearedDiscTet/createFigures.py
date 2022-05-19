import os
import subprocess
import casefoam
from typing import Union, List

curPath = os.getcwd()
os.makedirs("figures",exist_ok=True)
figPath = os.path.join(curPath,"figures")

for c in casefoam.of_cases("Cases"):
    os.chdir(c)
    l = c.split("/")
    c_name = "-".join(l[1:])
    image_name = f"vortexShearedDiscTet-{c_name}.png"
    image_path = os.path.join(curPath,"figures",image_name)
    subprocess.call(['pvbatch', 'stateFile_to_Image.py','interface.pvsm',image_path])
    os.chdir(curPath)
