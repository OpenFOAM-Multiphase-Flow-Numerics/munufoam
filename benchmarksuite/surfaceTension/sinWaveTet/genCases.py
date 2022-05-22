import casefoam

case = [['geoVoF_Brackbill','algVoF_Brackbill'],
        ['coarse', 'mid', 'fine']]

def update_Method(advect,recon,surfTenModel):
        return { 'system/simulationParameters': {'ADVECTSCHEME': advect,
                                                 'RECONSCHEME': recon,
                                                 'SURFTENMODEL': surfTenModel }}

def changeGeoFile(Nx):
    return {
                'triSquare.geo': 
                { '#!stringManipulation': 
                    {
                        "nx=20;": f"nx={Nx};",
                        "nz=60;": f"nz={3*Nx};"
                    }
                }
            }


data = {'geoVoF_Brackbill': update_Method("isoAdvection","plicRDF","Brackbill"),
        'algVoF_Brackbill': update_Method("MULESScheme","isoSurface","Brackbill"),
        'coarse': changeGeoFile(32),
        'mid': changeGeoFile(64),
        'fine': changeGeoFile(128)}

casefoam.mkCases('sinWaveTet', case, data, 'tree',writeDir='Cases')
