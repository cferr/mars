# Jacobi 2D, skewed rectangular tiling

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal

############################# Definitions ##############################

# Only parameter we'll want to change here
jacobi2d_r_TileSizes = [4, 5, 7]

jacobi2d_r_Space = ["t", "i", "j"]
jacobi2d_r_Dependences = [[1, 0, 1], [1, 1, 0], [1, 0, 0], [1, -1, 0], [1, 0, -1]]
jacobi2d_r_Hyperplanes = ["t", "t+i", "t+j"]
jacobi2d_r_ParamNames = ["k1", "k2", "k3"]
jacobi2d_r_NormalVectors = [[1, 0, 0], [1, 1, 0], [1, 0, 1]]

########################## Consistency checks ##########################

assert checkTilingLegal(
    jacobi2d_r_Dependences, jacobi2d_r_NormalVectors
), "Illegal tiling"

############################# Compute MARS #############################

jacobi2d_r_mars = getMARS(
    jacobi2d_r_Space,
    jacobi2d_r_Dependences,
    jacobi2d_r_Hyperplanes,
    jacobi2d_r_NormalVectors,
    jacobi2d_r_TileSizes,
)

jacobi2d_r_flowin_mars = getFlowInMARS(
    jacobi2d_r_Space,
    jacobi2d_r_Hyperplanes,
    jacobi2d_r_ParamNames,
    jacobi2d_r_TileSizes,
    jacobi2d_r_mars,
)
jacobi2d_r_flowout_mars = getFlowOutMARS(
    jacobi2d_r_Space,
    jacobi2d_r_Hyperplanes,
    jacobi2d_r_ParamNames,
    jacobi2d_r_TileSizes,
    jacobi2d_r_mars,
)
print(f"There are {len(jacobi2d_r_flowin_mars)} flow-in MARS, {len(jacobi2d_r_flowout_mars)} flow-out MARS per tile")
