# Canonical Deps 3D

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal

############################# Definitions ##############################

# Only parameter we'll want to change here
canonical3d_TileSizes = [10, 10, 10]

canonical3d_Space = ["k", "i", "j"]
canonical3d_Dependences = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
canonical3d_Hyperplanes = ["k", "i", "j"]
canonical3d_ParamNames = ["k1", "k2", "k3"]
canonical3d_NormalVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

########################## Consistency checks ##########################

assert checkTilingLegal(
    canonical3d_Dependences, canonical3d_NormalVectors
), "Illegal tiling"

############################# Compute MARS #############################

canonical3d_mars = getMARS(
    canonical3d_Space,
    canonical3d_Dependences,
    canonical3d_Hyperplanes,
    canonical3d_NormalVectors,
    canonical3d_TileSizes,
)

canonical3d_flowin_mars = getFlowInMARS(
    canonical3d_Space,
    canonical3d_Hyperplanes,
    canonical3d_ParamNames,
    canonical3d_TileSizes,
    canonical3d_mars,
)
canonical3d_flowout_mars = getFlowOutMARS(
    canonical3d_Space,
    canonical3d_Hyperplanes,
    canonical3d_ParamNames,
    canonical3d_TileSizes,
    canonical3d_mars,
)

print(f"There are {len(canonical3d_flowin_mars)} flow-in MARS, {len(canonical3d_flowout_mars)} flow-out MARS per tile")

