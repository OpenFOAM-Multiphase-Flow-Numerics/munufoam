import casefoam

case = [["gravity", "staticGravity","staticGravityPLIC","gravityRecon","gravityRecon2"],
        #["staticGravityPLIC"],
       ["Grid1", "Grid2", "Grid3"]]


def changeGravModel(model):
    return {"system/simulationParameters": {"GRAVITYMODEL": model}}

def update_Method(advect,recon,surfTenModel):
        return { 'system/simulationParameters': {'ADVECTSCHEME': advect,
                                                 'RECONSCHEME': recon,
                                                 'SURFTENMODEL': surfTenModel }}

def changeBlockMesh(Nx):
    return { 'system/simulationParameters': {'nx': str(Nx),
                                             'ny': str(3*Nx),
                                             'nz': "1"}}

data = {"gravity": changeGravModel("gravity"),
        "staticGravity": changeGravModel("staticGravity"),
        "staticGravityPLIC": changeGravModel("staticGravityPLIC"),
        "gravityRecon": changeGravModel("gravityRecon"),
        "gravityRecon2": changeGravModel("gravityRecon2"),
        'Grid1': changeBlockMesh(32),
        'Grid2': changeBlockMesh(64),
        'Grid3': changeBlockMesh(128)}

casefoam.mkCases('sinWave', case, data, 'tree',writeDir='Cases')
