import casefoam

case = [['isoAlpha', 'plicRDF','gradAlpha'],
        ['x8','x16','x32', 'x64', 'x128', 'x256', 'x512', 'x1024']]


def updateScheme(recon):
        return {
                    'system/simulationParameter': {"RECONSCHEME": recon }
               }


def changeBlockMesh(Nx):
    return {
                'system/simulationParameter': {
                    "nx": str(Nx),
                    "ny": "1",
                    "nz": str(Nx)
                }
            }

data = {'isoAlpha': updateScheme("isoAlpha"),
        'plicRDF': updateScheme("plicRDF"),
        'gradAlpha': updateScheme("gradAlpha"),
        'x8': changeBlockMesh(8),
        'x16': changeBlockMesh(16),
        'x32': changeBlockMesh(32),
        'x64': changeBlockMesh(64),
        'x128': changeBlockMesh(128),
        'x256': changeBlockMesh(256),
        'x512': changeBlockMesh(512),
        'x1024': changeBlockMesh(1024)}

casefoam.mkCases('circle', case, data, 'tree', writeDir='Cases')
