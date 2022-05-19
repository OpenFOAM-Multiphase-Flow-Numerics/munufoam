import casefoam

case = [['isoAlpha', 'plicRDF','gradAlpha'],
        ['x8','x16','x32', 'x64']]#, 'x128']] # x128 needs more than 16 Gb


def updateScheme(recon):
        return {
                    'system/simulationParameter': {"RECONSCHEME": recon }
               }


def changeGeoFile(Nx):
    return {
                'triBlock.geo': 
                { '#!stringManipulation': 
                    {
                        "nx=10;": f"nx={Nx};"
                    }
                }
            }

data = {'isoAlpha': updateScheme("isoAlpha"),
        'plicRDF': updateScheme("plicRDF"),
        'gradAlpha': updateScheme("gradAlpha"),
        'x8': changeGeoFile(8),
        'x16': changeGeoFile(16),
        'x32': changeGeoFile(32),
        'x64': changeGeoFile(64),
        'x128': changeGeoFile(128)}

casefoam.mkCases('sphereTet', case, data, 'tree', writeDir='Cases')
