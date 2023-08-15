import casefoam

case = [["gravity", "staticGravity","staticGravityPLIC"],
       ["Grid1", "Grid2", "Grid3"]]


def changeGravModel(model):
    return {"system/simulationParameters": {"GRAVITYMODEL": model}}

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

data = {"gravity": changeGravModel("gravity"),
        "staticGravity": changeGravModel("staticGravity"),
        "staticGravityPLIC": changeGravModel("staticGravityPLIC"),
        'Grid1': changeGeoFile(32),
        'Grid2': changeGeoFile(64),
        'Grid3': changeGeoFile(128)}

casefoam.mkCases('sinWaveTet', case, data, 'tree',writeDir='Cases')
