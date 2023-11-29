# GEMM

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal

############################# Definitions ##############################

# Only parameter we'll want to change here
gemm_TileSizes = [10, 20, 20]

gemm_Space = ["i", "k", "j"]
gemm_Dependences = [[0, 1, 0]]
gemm_Hyperplanes = ["i", "j", "k"]
gemm_ParamNames = ["k1", "k2", "k3"]
gemm_NormalVectors = [[1, 0, 0], [0, 0, 1], [0, 1, 0]]

########################## Consistency checks ##########################

assert checkTilingLegal(gemm_Dependences, gemm_NormalVectors), "Illegal tiling"

############################# Compute MARS #############################

gemm_mars = getMARS(
    gemm_Space, gemm_Dependences, gemm_Hyperplanes, gemm_NormalVectors, gemm_TileSizes
)

gemm_flowin_mars = getFlowInMARS(
    gemm_Space,
    gemm_Hyperplanes,
    gemm_ParamNames,
    gemm_TileSizes,
    gemm_mars,
)
gemm_flowout_mars = getFlowOutMARS(
    gemm_Space,
    gemm_Hyperplanes,
    gemm_ParamNames,
    gemm_TileSizes,
    gemm_mars,
)
print(f"There are {len(gemm_flowin_mars)} flow-in MARS, {len(gemm_flowout_mars)} flow-out MARS per tile")
