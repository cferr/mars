# Smith-Waterman

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal

############################# Definitions ##############################

# Only parameter we'll want to change here
sw_TileSizes = [100, 100]

sw_Space = ["i", "j"]
sw_Dependences = [[1, 0], [0, 1], [1, 1]]
sw_Hyperplanes = ["i+j", "j"]
sw_ParamNames = ["k1", "k2"]
sw_NormalVectors = [[1, 1], [0, 1]]

########################## Consistency checks ##########################

assert checkTilingLegal(sw_Dependences, sw_NormalVectors), "Illegal tiling"

############################# Compute MARS #############################

sw_mars = getMARS(
    sw_Space, sw_Dependences, sw_Hyperplanes, sw_NormalVectors, sw_TileSizes
)

sw_flowin_mars = getFlowInMARS(
    sw_Space, sw_Hyperplanes, sw_ParamNames, sw_TileSizes, sw_mars
)
sw_flowout_mars = getFlowOutMARS(
    sw_Space, sw_Hyperplanes, sw_ParamNames, sw_TileSizes, sw_mars
)
print(f"There are {len(sw_flowin_mars)} flow-in MARS, {len(sw_flowout_mars)} flow-out MARS per tile")

