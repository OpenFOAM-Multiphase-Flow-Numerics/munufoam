import casefoam
import numpy as np

case = [['geoVoF_Brackbill','algVoF_Brackbill'],
        ['Grid1' ,'Grid2','Grid3','Grid4','Grid5','Grid6','Grid7','Grid8','Grid9','Grid10']]

def update_Method(advect,recon,surfTenModel):
        return { 'system/simulationParameters': {'ADVECTSCHEME': advect,
                                                 'RECONSCHEME': recon,
                                                 'SURFTENMODEL': surfTenModel }}

def changeGeoFile(Nx):
    return {
                'triSquare.geo': 
                { '#!stringManipulation': 
                    {
                        "nx=40;": f"nx={4*Nx};",
                        "nz=10;": f"nz={Nx};"
                    }
                }
            }

res = np.linspace(10,128,10)
res = res.astype(int)

data = {'geoVoF_Brackbill': update_Method("isoAdvection","plicRDF","Brackbill"),
        'algVoF_Brackbill': update_Method("MULESScheme","isoSurface","Brackbill"),
        'Grid1': changeGeoFile(res[0]),
        'Grid2': changeGeoFile(res[1]),
        'Grid3': changeGeoFile(res[2]),
        'Grid4': changeGeoFile(res[3]),
        'Grid5': changeGeoFile(res[4]),
        'Grid6': changeGeoFile(res[5]),
        'Grid7': changeGeoFile(res[6]),
        'Grid8': changeGeoFile(res[7]),
        'Grid9': changeGeoFile(res[8]),
        'Grid10': changeGeoFile(res[9])}

casefoam.mkCases('advectedCircleTet', case, data, 'tree',writeDir='Cases')
