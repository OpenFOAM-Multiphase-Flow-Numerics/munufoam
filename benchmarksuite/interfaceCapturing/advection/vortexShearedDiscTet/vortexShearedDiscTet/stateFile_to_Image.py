from paraview.simple import *
import os, sys

if len(sys.argv) < 2:
    sys.exit('Two parameters, name of an input file and the output path')
if not os.path.exists(sys.argv[1]):
    sys.exit('Specified file "' + sys.argv[1] + '" does not exists!')

stateFile = sys.argv[1]
image_file = sys.argv[2]


servermanager.LoadState(stateFile)
SetActiveView(GetRenderView())
view = GetActiveView()

scene = GetAnimationScene()

print("image file path: ",image_file)

# save animation images
SaveScreenshot(image_file,view)

