# Jacobi 2D, diamond tiling

from mars.compute.mars import getMARS, getFlowInMARS, getFlowOutMARS
from mars.polyhedral.utils import checkTilingLegal
from islpy import Set

############################# Definitions ##############################

jacobi2d_d_Space = ["t", "i", "j"]
jacobi2d_d_Dependences = [[1, 0, 1], [1, 1, 0], [1, 0, 0], [1, -1, 0], [1, 0, -1]]
jacobi2d_d_Hyperplanes = ["t+i", "t+j", "t-i", "t-j"]
jacobi2d_d_ParamNames = ["k1", "k2", "k3", "k4"]
jacobi2d_d_NormalVectors = [[1, 1, 0], [1, 0, 1], [1, -1, 0], [1, 0, -1]]
jacobi2d_d_TileSizes = [20, 20, 20, 20]

# I do not know how to automatically infer these... but let's say there is a way.
jacobi2d_d_ParamRelations = [
    Set("[k1, k2, k3, k4] -> { : k4 = k1 - k2 + k3 }"),
    Set("[k1, k2, k3, k4] -> { : k4 = -1 + k1 - k2 + k3 }"),
    Set("[k1, k2, k3, k4] -> { : k4 = 1 + k1 - k2 + k3 }"),
]

jacobi2d_d_ParamRelationNames = {
    jacobi2d_d_ParamRelations[0]: "R0",
    jacobi2d_d_ParamRelations[1]: "R1",
    jacobi2d_d_ParamRelations[2]: "R2",
}

########################## Consistency checks ##########################

assert checkTilingLegal(
    jacobi2d_d_Dependences, jacobi2d_d_NormalVectors
), "Illegal tiling"

############################# Compute MARS #############################

# Compute MARS
jacobi2d_d_mars = getMARS(
    jacobi2d_d_Space,
    jacobi2d_d_Dependences,
    jacobi2d_d_Hyperplanes,
    jacobi2d_d_NormalVectors,
    jacobi2d_d_TileSizes,
)


jacobi2d_d_flowin_mars = {}
jacobi2d_d_flowout_mars = {}

for relation in jacobi2d_d_ParamRelations:
    # Gist MARS with context
    jacobi2d_d_flowin_mars[relation] = getFlowInMARS(
        jacobi2d_d_Space,
        jacobi2d_d_Hyperplanes,
        jacobi2d_d_ParamNames,
        jacobi2d_d_TileSizes,
        jacobi2d_d_mars,
        ParamRelation=relation,
    )
    jacobi2d_d_flowout_mars[relation] = getFlowOutMARS(
        jacobi2d_d_Space,
        jacobi2d_d_Hyperplanes,
        jacobi2d_d_ParamNames,
        jacobi2d_d_TileSizes,
        jacobi2d_d_mars,
        ParamRelation=relation,
    )
    print(f"There are {len(jacobi2d_d_flowin_mars[relation])} flow-in MARS, {len(jacobi2d_d_flowout_mars[relation])} flow-out MARS per tile with relation {relation}")

