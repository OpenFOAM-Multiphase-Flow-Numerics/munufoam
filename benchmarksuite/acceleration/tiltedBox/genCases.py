import casefoam

# directory of the base case
baseCase = "tiltedBox"

# list of parent, child and grandchild names
caseStructure = [#["gravity", "staticGravity"],
                 ["gravity", "staticGravity","staticGravityPLIC"],
                 ["Grid1", "Grid2", "Grid3", "Grid4"]]


def changeGravModel(model):
    return {"system/simulationParameters": {"GRAVITYMODEL": model}}


# this function does the same as update_coarse etc but is more elegant
def changeResolution(Nx):
    return {
        "system/simulationParameters": {"nx": str(Nx),"nz":  str(Nx)}
    }


# dictionary of data to update
caseData = {
    "gravity": changeGravModel("gravity"),
    "staticGravity": changeGravModel("staticGravity"),
    "staticGravityPLIC": changeGravModel("staticGravityPLIC"),
    "Grid1": changeResolution(32),
    "Grid2": changeResolution(64),
    "Grid3": changeResolution(128),
    "Grid4": changeResolution(256),
    "Grid5": changeResolution(512),
}

# generate cases
casefoam.mkCases(baseCase, caseStructure, caseData, hierarchy="tree", writeDir="Cases")
