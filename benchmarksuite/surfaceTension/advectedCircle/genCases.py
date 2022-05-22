import casefoam
import numpy as np

case = [['geoVoF_Brackbill','algVoF_Brackbill'],
        ['Grid1' ,'Grid2','Grid3','Grid4','Grid5','Grid6','Grid7','Grid8','Grid9','Grid10']]

def update_Method(advect,recon,surfTenModel):
        return { 'system/simulationParameters': {'ADVECTSCHEME': advect,
                                                 'RECONSCHEME': recon,
                                                 'SURFTENMODEL': surfTenModel }}

def changeBlockMesh(Nx):
    return { 'system/simulationParameters': {'nx': str(4*Nx),
                                             'ny': "1",
                                             'nz': str(Nx)}}

res = np.linspace(10,128,10)
res = res.astype(int)

data = {'geoVoF_Brackbill': update_Method("isoAdvection","plicRDF","Brackbill"),
        'algVoF_Brackbill': update_Method("MULESScheme","isoSurface","Brackbill"),
        'Grid1': changeBlockMesh(res[0]),
        'Grid2': changeBlockMesh(res[1]),
        'Grid3': changeBlockMesh(res[2]),
        'Grid4': changeBlockMesh(res[3]),
        'Grid5': changeBlockMesh(res[4]),
        'Grid6': changeBlockMesh(res[5]),
        'Grid7': changeBlockMesh(res[6]),
        'Grid8': changeBlockMesh(res[7]),
        'Grid9': changeBlockMesh(res[8]),
        'Grid10': changeBlockMesh(res[9])}

casefoam.mkCases('advectedCircle', case, data, 'tree',writeDir='Cases')
