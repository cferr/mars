# Jacobi 1D

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal

############################# Definitions ##############################

# Only parameter we'll want to change here
jacobi1d_TileSizes = [45000, 45000]

jacobi1d_Space = ["t", "i"]
jacobi1d_Dependences = [[1, -1], [1, 0], [1, 1]]
jacobi1d_Hyperplanes = ["t + i", "t - i"]
jacobi1d_ParamNames = ["k1", "k2"]
jacobi1d_NormalVectors = [[1, 1], [1, -1]]

########################## Consistency checks ##########################

# Check that tiling is legal (more restrictive than actually needed, think of height-1 tiles)
assert checkTilingLegal(jacobi1d_Dependences, jacobi1d_NormalVectors), "Illegal tiling"

############################# Compute MARS #############################

jacobi1d_mars = getMARS(
    jacobi1d_Space,
    jacobi1d_Dependences,
    jacobi1d_Hyperplanes,
    jacobi1d_NormalVectors,
    jacobi1d_TileSizes,
)
jacobi1d_flowin_mars = getFlowInMARS(
    jacobi1d_Space,
    jacobi1d_Hyperplanes,
    jacobi1d_ParamNames,
    jacobi1d_TileSizes,
    jacobi1d_mars,
)
jacobi1d_flowout_mars = getFlowOutMARS(
    jacobi1d_Space,
    jacobi1d_Hyperplanes,
    jacobi1d_ParamNames,
    jacobi1d_TileSizes,
    jacobi1d_mars,
)
print(f"There are {len(jacobi1d_flowin_mars)} flow-in MARS, {len(jacobi1d_flowout_mars)} flow-out MARS per tile")

