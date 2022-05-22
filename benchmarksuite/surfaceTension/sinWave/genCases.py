import casefoam

case = [['geoVoF_Brackbill','algVoF_Brackbill'],
        ['coarse', 'mid', 'fine']]

def update_Method(advect,recon,surfTenModel):
        return { 'system/simulationParameters': {'ADVECTSCHEME': advect,
                                                 'RECONSCHEME': recon,
                                                 'SURFTENMODEL': surfTenModel }}

def changeBlockMesh(Nx):
    return { 'system/simulationParameters': {'nx': str(Nx),
                                             'ny': str(3*Nx),
                                             'nz': "1"}}

data = {'geoVoF_Brackbill': update_Method("isoAdvection","plicRDF","Brackbill"),
        'algVoF_Brackbill': update_Method("MULESScheme","isoSurface","Brackbill"),
        'coarse': changeBlockMesh(32),
        'mid': changeBlockMesh(64),
        'fine': changeBlockMesh(128)}

casefoam.mkCases('sinWave', case, data, 'tree',writeDir='Cases')
