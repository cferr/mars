# Seidel 2D

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal

############################# Definitions ##############################

# Only parameter we'll want to change here
seidel2d_TileSizes = [4, 10, 10]

seidel2d_Space = ["t", "i", "j"]
seidel2d_Dependences = [
    [0, 1, 1],
    [0, 0, 1],
    [1, -1, 1],
    [0, 1, 0],
    [1, 0, 0],
    [1, -1, 0],
    [0, 1, -1],
    [1, 0, -1],
    [1, -1, -1],
]
seidel2d_Hyperplanes = ["t", "t + i", "4t + 2i + j"]
seidel2d_ParamNames = ["k1", "k2", "k3"]
seidel2d_NormalVectors = [[1, 0, 0], [1, 1, 0], [4, 2, 1]]

########################## Consistency checks ##########################

assert checkTilingLegal(seidel2d_Dependences, seidel2d_NormalVectors), "Illegal tiling"

############################# Compute MARS #############################

seidel2d_mars = getMARS(
    seidel2d_Space,
    seidel2d_Dependences,
    seidel2d_Hyperplanes,
    seidel2d_NormalVectors,
    seidel2d_TileSizes,
)

seidel2d_flowin_mars = getFlowInMARS(
    seidel2d_Space,
    seidel2d_Hyperplanes,
    seidel2d_ParamNames,
    seidel2d_TileSizes,
    seidel2d_mars,
)

seidel2d_flowout_mars = getFlowOutMARS(
    seidel2d_Space,
    seidel2d_Hyperplanes,
    seidel2d_ParamNames,
    seidel2d_TileSizes,
    seidel2d_mars,
)

print(f"There are {len(seidel2d_flowin_mars)} flow-in MARS, {len(seidel2d_flowout_mars)} flow-out MARS per tile")

