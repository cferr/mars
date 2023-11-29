from islpy import AstBuild, Map, Set, UnionMap, dim_type
from mars.compute.mars import getConsumerTiles
from mars.polyhedral.utils import getTileInnerSet, getDependencesAsUnionMap
from mars.polyhedral.barvinok import getBarvinokCardAsPWQPolynomial

from functools import reduce
import re


def genRelationComments(Relations):
    comments = []
    for R in Relations:
        rel_idx = Relations.index(R)
        comments += [f"// Relation R{rel_idx} is : {R.to_str()}"]
    return "\n".join(comments)


def getPolyhedralCodegenMacros():
    return (
        "#define ceild(n,d) (((n) < 0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))\n"
        "#define floord(x,y) (((x) < 0)? -((-(x)+(y)-1)/(y)) : (x)/(y))\n"
        "#define max(x,y)    ((x) > (y)? (x) : (y))\n"
        "#define min(x,y)    ((x) < (y)? (x) : (y))\n"
    )


def wrapIntoFunction(FunctionType, FunctionName, FunctionArgs, FunctionBody):
    return (
        str(
            FunctionType
            + " "
            + FunctionName
            + "("
            + ", ".join(
                [ArgType + " " + ArgName for (ArgType, ArgName) in FunctionArgs]
            )
            + ") {\n"
            + FunctionBody
        ).replace("\n", "\n  ")
        + "\n}"
    )


def generateTileSkeleton(
    Space,
    Dependences,
    Hyperplanes,
    ParamNames,
    TileSizes,
    StmtToSpaceMap=None,
    StmtRemap=None,
    StmtNamePrepend=None,
    IntraTileSchedule=None,
    all_rel=None,
):
    Decal = [0 for i in range(len(Hyperplanes))]
    TileSpace = getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal)
    if not all_rel is None:
        TileSpace = TileSpace.gist_params(all_rel)

    if not StmtRemap is None:
        TileSpace = TileSpace.apply(StmtRemap.copy())

    if IntraTileSchedule is None:
        # TileSchedule = TileSpace.identity() # STOP! This schedule may be plain ILLEGAL!
        TileSchedule = TileSpace.compute_schedule(
            getDependencesAsUnionMap(Space, Dependences), UnionMap("{ }")
        )
        TileScheduleMap = TileSchedule.get_map().intersect_domain(
            TileSchedule.get_domain()
        )
    else:
        TileScheduleMap = IntraTileSchedule.intersect_domain(TileSpace)

    if StmtToSpaceMap is None:
        TileScheduleStmt = TileScheduleMap
    else:
        TileScheduleStmt = TileScheduleMap.apply_domain(StmtToSpaceMap.copy().reverse())

    # Rename statements here if necessary (e.g., to add the "HOST_" that triggers boundary checks)
    if not StmtNamePrepend is None:
        prependedMaps = []
        TileScheduleStmtMaps = TileScheduleStmt.get_map_list()
        for i in range(TileScheduleStmtMaps.n_map()):
            TileScheduleStmtMap = TileScheduleStmtMaps.get_map(i)
            OriginalName = TileScheduleStmtMap.get_tuple_name(dim_type.in_)
            if OriginalName is None:
                OriginalName = ""
            RenamedStmtMap = TileScheduleStmtMap.set_tuple_name(
                dim_type.in_, StmtNamePrepend + OriginalName
            )
            prependedMaps += [RenamedStmtMap]
        # Make a union map back from the maps with prepended names
        TileScheduleStmt = reduce(
            lambda a, b: a.union(b), [x.to_union_map() for x in prependedMaps]
        )

    Ctx = TileScheduleStmt.get_ctx()
    islAstBuild = AstBuild.alloc(Ctx)
    AST = islAstBuild.ast_from_schedule(TileScheduleStmt)
    coreCode = "  " + AST.to_C_str().replace("\n", "\n  ")[:-2]

    functionArgs = [("int", P) for P in ParamNames]

    #     functionMockup = "void tile_" + Name + "(" + ", ".join([f"int {P}" for P in ParamNames]) +") {\n" + \
    #         "  // You should use a scheduler like PLuTo to get parallelism from this code.\n" + \
    #         coreCode + "\n}"
    #
    #     return functionMockup
    return (coreCode, functionArgs)


def generateTileSkeletonWithRelations(
    Name,
    Space,
    Dependences,
    Hyperplanes,
    ParamNames,
    TileSizes,
    Relations,
    StmtToSpaceMap=None,
    StmtRemap=None,
):
    Decal = [0 for i in range(len(Hyperplanes))]
    TileSpace = getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal)
    if not StmtRemap is None:
        TileSpace = TileSpace.apply(StmtRemap.copy())

    functions = []
    for R in Relations:
        rel_idx = Relations.index(R)
        # TileSchedule = TileSpace.identity() # STOP! This schedule may be plain ILLEGAL!
        TileSchedule = (
            TileSpace.copy()
            .gist_params(R.copy())
            .compute_schedule(
                getDependencesAsUnionMap(Space, Dependences), UnionMap("{ }")
            )
        )
        TileScheduleMap = TileSchedule.get_map().intersect_domain(
            TileSchedule.get_domain()
        )

        if StmtToSpaceMap is None:
            TileScheduleStmt = TileScheduleMap
        else:
            TileScheduleStmt = TileScheduleMap.apply_domain(
                StmtToSpaceMap.copy().reverse()
            )

        islAstBuild = AstBuild.from_context(R)
        AST = islAstBuild.ast_from_schedule(TileScheduleStmt)
        coreCode = "  " + AST.to_C_str().replace("\n", "\n  ")[:-2]

        functionMockup = (
            "void tile_"
            + Name
            + f"_R{rel_idx}("
            + ", ".join([f"int {P}" for P in ParamNames])
            + ") {\n"
            + "  // You should use a scheduler like PLuTo to get parallelism from this code.\n"
            + coreCode
            + "\n}"
        )
        functions += [functionMockup]

    return "\n".join(functions)


def genRelationDispatchCode(Relations, RelationNames=[]):
    if len(RelationNames) > 0:
        assert len(RelationNames) == len(Relations), "One name per relation"
    sp = Relations[0].get_space()
    UM = UnionMap.empty_space(sp)
    for R in Relations:
        rel_idx = Relations.index(R)
        name_suffix = (
            RelationNames[rel_idx] if len(RelationNames) > 0 else f"R{rel_idx}"
        )
        Nparams = R.dim(dim_type.param)
        Rset = Map.from_domain(
            R.copy()
            .move_dims(dim_type.set, 0, dim_type.param, 0, Nparams)
            .identity()
            .move_dims(dim_type.param, 0, dim_type.in_, 0, Nparams)
            .range()
        ).set_tuple_name(dim_type.in_, name_suffix)
        UM = UM.add_map(Rset)
    islASTBuild = AstBuild.alloc(UM.get_ctx())
    ast = islASTBuild.ast_from_schedule(UM)
    dispatchCode = "  " + ast.to_C_str().replace("\n", "\n  ")[:-2]
    return dispatchCode


def genRelationDispatchFunction(TargetFunctionName, Relations, RelationNames=[]):
    dispatchCode = genRelationDispatchCode(Relations, RelationNames)
    sp = Relations[0].get_space()
    params = [
        "int " + sp.get_dim_name(dim_type.param, i)
        for i in range(sp.dim(dim_type.param))
    ]
    return (
        "void "
        + TargetFunctionName
        + "("
        + ", ".join(params)
        + ") {\n"
        + dispatchCode
        + "}"
    )


def getAllTiles(Space, Hyperplanes, ParamNames, TileSizes, Boundaries):
    # Construct hyperplane boundaries (inequalities)
    BoundariesAsISL = []
    for P, H, E in Boundaries:
        SetStr = ""
        if len(P) > 0:
            SetStr += "[" + ", ".join(P) + "] -> "
        NormalEq = (
            "("
            + " + ".join(
                [str(H[i]) + " * " + Space[i] for i in range(len(H)) if H[i] != 0]
            )
            + ")"
        )
        SetStr += "{ [" + ", ".join(Space) + "] : " + NormalEq + " >= " + E + " }"
        BoundariesAsISL += [Set(SetStr)]
    # Get a tile
    Decal = [0 for _ in range(len(Hyperplanes))]
    Tile = getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal)
    AllBoundariesIntersect = Tile.copy()
    for B in BoundariesAsISL:
        AllBoundariesIntersect = AllBoundariesIntersect.intersect(B).coalesce()
    Result = AllBoundariesIntersect.project_out(
        dim_type.set, 0, AllBoundariesIntersect.dim(dim_type.set)
    ).move_dims(
        dim_type.set,
        0,
        dim_type.param,
        AllBoundariesIntersect.dim(dim_type.param) - len(Hyperplanes),
        len(Hyperplanes),
    )
    return Result


# All tiles that do not contain a border of the space.
def getInset(Space, Hyperplanes, ParamNames, TileSizes, Boundaries):
    # Construct hyperplane boundaries (equalities)
    BoundariesAsISL = []
    for P, H, E in Boundaries:
        SetStr = ""
        if len(P) > 0:
            SetStr += "[" + ", ".join(P) + "] -> "
        NormalEq = (
            "("
            + " + ".join(
                [str(H[i]) + " * " + Space[i] for i in range(len(H)) if H[i] != 0]
            )
            + ")"
        )
        SetStr += "{ [" + ", ".join(Space) + "] : " + NormalEq + " = " + E + " }"
        BoundariesAsISL += [Set(SetStr)]
    # Get a tile
    Decal = [0 for _ in range(len(Hyperplanes))]
    Tile = getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal)
    AllBoundariesUnion = Set.empty(Tile.get_space())
    for B in BoundariesAsISL:
        AllBoundariesUnion = AllBoundariesUnion.union(B).coalesce()
    FiniteBorders = (
        Tile.intersect(AllBoundariesUnion)
        .project_out(dim_type.set, 0, AllBoundariesUnion.dim(dim_type.set))
        .move_dims(
            dim_type.set,
            0,
            dim_type.param,
            AllBoundariesUnion.dim(dim_type.param) - len(Hyperplanes),
            len(Hyperplanes),
        )
    )
    return getAllTiles(Space, Hyperplanes, ParamNames, TileSizes, Boundaries).subtract(
        FiniteBorders
    )


def generateInsetTestingCode(Space, Hyperplanes, ParamNames, TileSizes, Boundaries):
    inset = getInset(Space, Hyperplanes, ParamNames, TileSizes, Boundaries)
    insetParam = (
        inset.move_dims(
            dim_type.param, inset.dim(dim_type.param), dim_type.set, 0, len(Hyperplanes)
        )
        .set_tuple_name("return_1")
        .identity()
    )
    islASTBuild = AstBuild.alloc(insetParam.get_ctx())
    ast = islASTBuild.ast_from_schedule(UnionMap.from_map(insetParam))
    cCode = ast.to_C_str()
    testProgStr = ""
    searchRes = re.search("if[ ]*\((.*)\)", cCode)
    if not searchRes is None:
        testProgStr = searchRes.group(1)
    return testProgStr


def generateInsetTestingFunction(Name, Space, Hyperplanes, TileSizes, Boundaries):
    inset = getInset(Space, Hyperplanes, TileSizes, Boundaries)
    insetParam = (
        inset.move_dims(
            dim_type.param, inset.dim(dim_type.param), dim_type.set, 0, len(Hyperplanes)
        )
        .set_tuple_name("return_1")
        .identity()
    )
    islASTBuild = AstBuild.alloc(insetParam.get_ctx())
    ast = islASTBuild.ast_from_schedule(UnionMap.from_map(insetParam))
    testProgStr = (
        ast.to_C_str().replace("\n", "\n  ")[:-2].replace("return_1()", "return 1")
        + "\n  return 0;\n"
    )

    params = [
        insetParam.get_dim_name(dim_type.param, i)
        for i in range(insetParam.dim(dim_type.param))
    ]
    arguments = [("int", p) for p in params]
    return wrapIntoFunction("int", Name, arguments, testProgStr)


# All tiles that contain some border of the space.
def getOutset(Space, Hyperplanes, ParamNames, TileSizes, Boundaries):
    # Construct hyperplane boundaries (equalities)
    BoundariesAsISL = []
    for P, H, E in Boundaries:
        SetStr = ""
        if len(P) > 0:
            SetStr += "[" + ", ".join(P) + "] -> "
        NormalEq = (
            "("
            + " + ".join(
                [str(H[i]) + " * " + Space[i] for i in range(len(H)) if H[i] != 0]
            )
            + ")"
        )
        SetStr += "{ [" + ", ".join(Space) + "] : " + NormalEq + " = " + E + " }"
        BoundariesAsISL += [Set(SetStr)]
    # Get a tile
    Decal = [0 for _ in range(len(Hyperplanes))]
    Tile = getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal)
    AllBoundariesUnion = Set.empty(Tile.get_space())
    for B in BoundariesAsISL:
        AllBoundariesUnion = AllBoundariesUnion.union(B).coalesce()
    FiniteBorders = (
        Tile.intersect(AllBoundariesUnion)
        .project_out(dim_type.set, 0, AllBoundariesUnion.dim(dim_type.set))
        .move_dims(
            dim_type.set,
            0,
            dim_type.param,
            AllBoundariesUnion.dim(dim_type.param) - len(Hyperplanes),
            len(Hyperplanes),
        )
    )
    return getAllTiles(Space, Hyperplanes, ParamNames, TileSizes, Boundaries).intersect(
        FiniteBorders
    )


def getFlowInPlusTile(Space, Dependences, Hyperplanes, ParamNames, TileSizes):
    currentTile = getTileInnerSet(
        Space,
        Hyperplanes,
        ParamNames,
        TileSizes,
        [0 for i in Hyperplanes],
    )
    SpaceStr = "[" + ", ".join(Space) + "]"
    # Make ISL map with dependences (reversed)
    DepMap = Map.empty(currentTile.copy().identity().get_space())
    for dep in Dependences:
        MapTarget = (
            "["
            + ", ".join(
                [
                    Space[i] + (" - " if dep[i] > 0 else " + ") + str(abs(dep[i]))
                    for i in range(len(dep))
                ]
            )
            + "]"
        )
        thisDepMap = Map("{ " + SpaceStr + " -> " + MapTarget + " }")
        DepMap = DepMap.union(thisDepMap)
    # Compute the flow-in + the tile (union instead of subtract)
    flowInParam = currentTile.copy().apply(DepMap).union(currentTile)
    return flowInParam


def getISLIterationSpace(Space, Boundaries):
    # Make a sample ISL space
    sampleISL = Set("{ [" + ", ".join(Space) + "] }")
    # Get ISL boundaries
    BoundariesAsISL = []
    for P, H, E in Boundaries:
        SetStr = ""
        if len(P) > 0:
            SetStr += "[" + ", ".join(P) + "] -> "
        NormalEq = (
            "("
            + " + ".join(
                [str(H[i]) + " * " + Space[i] for i in range(len(H)) if H[i] != 0]
            )
            + ")"
        )
        SetStr += "{ [" + ", ".join(Space) + "] : " + NormalEq + " >= " + E + " }"
        BoundariesAsISL += [Set(SetStr)]
    AllBoundariesIntersection = Set.universe(sampleISL.get_space())
    # Intersect flow-in with boundaries, and get corresponding tiles
    for B in BoundariesAsISL:
        AllBoundariesIntersection = AllBoundariesIntersection.intersect(B).coalesce()
    return AllBoundariesIntersection


# Test if a tile's flow-in does not cross the boundaries
def getFPGAValidDomain(
    Space, Dependences, Hyperplanes, ParamNames, TileSizes, Boundaries
):
    currentTile = getTileInnerSet(
        Space,
        Hyperplanes,
        ParamNames,
        TileSizes,
        [0 for i in Hyperplanes],
    )
    flowInParam = getFlowInPlusTile(
        Space, Dependences, Hyperplanes, ParamNames, TileSizes
    )
    # Get ISL boundaries
    BoundariesAsISL = []
    for P, H, E in Boundaries:
        SetStr = ""
        if len(P) > 0:
            SetStr += "[" + ", ".join(P) + "] -> "
        NormalEq = (
            "("
            + " + ".join(
                [str(H[i]) + " * " + Space[i] for i in range(len(H)) if H[i] != 0]
            )
            + ")"
        )
        SetStr += "{ [" + ", ".join(Space) + "] : " + NormalEq + " = " + E + " }"
        BoundariesAsISL += [Set(SetStr)]
    AllBoundariesUnion = Set.empty(currentTile.get_space())
    # Intersect flow-in with boundaries, and get corresponding tiles
    for B in BoundariesAsISL:
        AllBoundariesUnion = AllBoundariesUnion.union(B).coalesce()
    FiniteBorders = (
        flowInParam.intersect(AllBoundariesUnion)
        .project_out(dim_type.set, 0, AllBoundariesUnion.dim(dim_type.set))
        .move_dims(
            dim_type.set,
            0,
            dim_type.param,
            AllBoundariesUnion.dim(dim_type.param) - len(Hyperplanes),
            len(Hyperplanes),
        )
    )
    return getAllTiles(Space, Hyperplanes, ParamNames, TileSizes, Boundaries).subtract(
        FiniteBorders
    )


# Get those tiles executed on FPGA whose consumers are all executed on FPGA
def getFPGAInnerDomain(
    Space, Dependences, Hyperplanes, NormalVectors, ParamNames, TileSizes, Boundaries
):
    FPGATiles = getFPGAValidDomain(
        Space, Dependences, Hyperplanes, ParamNames, TileSizes, Boundaries
    )
    innerDomainTiles = FPGATiles
    # Get consumer tiles
    consumers = getConsumerTiles(
        Space, Dependences, Hyperplanes, NormalVectors, TileSizes
    )
    for cons in consumers:
        # Get consumer as ISL
        interTileDep = [0 for _ in Hyperplanes]
        for x in cons:
            interTileDep[x] = 1
        islDeps = getDependencesAsUnionMap(
            [f"i{i}" for i in range(len(Hyperplanes))], [interTileDep]
        )
        subjectTiles = (
            islDeps.as_map()
            .intersect_domain(FPGATiles)
            .intersect_range(FPGATiles)
            .domain()
        )
        innerDomainTiles = innerDomainTiles.intersect(subjectTiles)
    return innerDomainTiles


# Generate code that browses all tiles for the space, distinguishing full ones (accelerated)
# from non-full ones (not accelerated).
def generateTileSchedule(
    Space,
    Hyperplanes,
    ParamNames,
    NormalVectors,
    TileSizes,
    Dependences,
    Boundaries,
    ParamContext=None,
):
    # Get all tiles, getInset, getOutset
    allTiles = getAllTiles(Space, Hyperplanes, ParamNames, TileSizes, Boundaries)
    # getInset = getInset(Space, Hyperplanes, TileSizes, Boundaries)
    # getOutset = getOutset(Space, Hyperplanes, TileSizes, Boundaries)
    # Parametrize getInset, getOutset
    # PInset = getInset.move_dims(dim_type.param, getInset.dim(dim_type.param), dim_type.set, 0, len(Hyperplanes))
    # POutset = getOutset.move_dims(dim_type.param, getOutset.dim(dim_type.param), dim_type.set, 0, len(Hyperplanes))
    # Generate dispatching function
    # dispatcher = genRelationDispatchFunction("tileDispatch_" + BenchmarkName, [PInset, POutset], \
    #                                 RelationNames = [BenchmarkName + "_Inset", BenchmarkName + "_Outset"])
    # Schedule all tiles according to dependences
    prodToConsDeps = getConsumerTiles(
        Space, Dependences, Hyperplanes, NormalVectors, TileSizes
    )  # form : each hyperplane which index is there is crossed

    # Convert the consumer tile hyperplane indices into dependence vectors
    interTileDeps = []
    for T in prodToConsDeps:
        DepVec = [0 for _ in range(len(Hyperplanes))]
        for H in T:
            DepVec[H] = 1
        interTileDeps += [DepVec]
    # Build an ISL representation of these dependence vectors
    interTileDepISL = UnionMap.empty(allTiles.get_space())
    for D in interTileDeps:
        addDep = Map(
            "{ ["
            + ", ".join(ParamNames)
            + "] -> ["
            + ", ".join(
                [p + (" + 1" if d == 1 else "") for (p, d) in zip(ParamNames, D)]
            )
            + "] }"
        )
        interTileDepISL = interTileDepISL.add_map(addDep)

    # print("Inter-tile dependence: " + str(interTileDepISL))

    # Schedule tiles according to interTileDep and build an AST
    interTileSchedule = allTiles.compute_schedule(
        interTileDepISL, UnionMap.empty(interTileDepISL.get_space())
    )
    interTileScheduleMap = interTileSchedule.get_map().intersect_domain(
        interTileSchedule.get_domain()
    )

    # Get the space size's parameters
    spaceParams = [
        allTiles.get_dim_name(dim_type.param, i)
        for i in range(allTiles.dim(dim_type.param))
    ]

    interTileScheduleMapNamed = UnionMap.empty_space(interTileScheduleMap.get_space())
    for i in range(interTileScheduleMap.n_map()):
        indivMap = interTileScheduleMap.get_map_list().get_at(i)
        addParamsMap = Map(
            "["
            + ", ".join(spaceParams)
            + "] -> { ["
            + ", ".join(ParamNames)
            + "] -> ["
            + ", ".join(spaceParams + ParamNames)
            + "] }"
        )
        indivMap = (
            indivMap.apply_domain(addParamsMap)
            .project_out(dim_type.in_, 0, len(spaceParams))
            .set_tuple_name(dim_type.in_, "TILE_DISPATCH_CALL")
        )
        interTileScheduleMapNamed = interTileScheduleMapNamed.add_map(indivMap)

    if ParamContext is None:
        ParamContextISL = Set("{ : }")
    else:
        ParamContextISL = ParamContext

    print(f"Tile schedule: {interTileScheduleMapNamed}")

    astBuild = AstBuild.from_context(ParamContextISL)
    tileBrowsingAST = (
        "  "
        + astBuild.ast_from_schedule(interTileScheduleMapNamed)
        .to_C_str()
        .replace("\n", "\n  ")[:-2]
    )

    # Finally wrap that into a function...
    #     tileBrowsingWrapperFun = "void browseTiles_" + BenchmarkName + "(" + ", ".join([f"int {P}" for P in spaceParams]) + ") {\n" + tileBrowsingAST + "}"

    #     return dispatcher + "\n" + tileBrowsingWrapperFun
    # TODO do params
    return (tileBrowsingAST, [])


def generateTileAllocationParams(Space, Hyperplanes, ParamNames, TileSizes, Boundaries):
    tileCoordSpace = getAllTiles(Space, Hyperplanes, ParamNames, TileSizes, Boundaries)
    NbTilesPerHyp = []
    MARSSpaceShiftPerHyp = []
    for i in range(len(TileSizes)):
        # Make a projection to only that dimension
        projMap = Map("{ [" + ", ".join(ParamNames) + f"] -> [{ParamNames[i]}]" + " }")
        allTilesOnThatAxis = tileCoordSpace.copy().apply(projMap)
        allTilesPWQ = getBarvinokCardAsPWQPolynomial(allTilesOnThatAxis)
        pwqPieces = allTilesPWQ.get_pieces()
        if len(pwqPieces) > 1:
            raise ValueError(
                "True Piecewise quasi-polynomial for cardinality of MARS space, unsupported"
            )
        (_, poly) = pwqPieces[0]
        polyStr = str(poly)
        searchRes = re.search(r"{ (.*) }", polyStr)
        if not searchRes is None:
            formula = searchRes.group(1)
            formulaMul = re.sub(r"([0-9]+)([a-zA-Z]+)", "\\1 * \\2", formula)
            NbTilesPerHyp += [formulaMul]
        else:
            raise ValueError("Failed to find polynomial")

        # Count how many tiles have a negative coordinate on that axis (we'll need to shift by that much)
        negativeSet = Set("{ [x] : x < 0 }")
        negativeSetCard = getBarvinokCardAsPWQPolynomial(
            allTilesOnThatAxis.intersect(negativeSet)
        )
        pwqPieces = negativeSetCard.get_pieces()
        if len(pwqPieces) > 1:
            raise ValueError(
                "True Piecewise quasi-polynomial for cardinality of MARS space, unsupported"
            )
        elif len(pwqPieces) == 0:
            # There are no negative tiles here
            MARSSpaceShiftPerHyp += ["0"]
        else:
            (_, poly) = pwqPieces[0]
            polyStr = str(poly)
            searchRes = re.search(r"{ (.*) }", polyStr)
            if not searchRes is None:
                formula = searchRes.group(1)
                formulaMul = re.sub(r"([0-9]+)([a-zA-Z]+)", "\\1 * \\2", formula)
                MARSSpaceShiftPerHyp += [formulaMul]
            else:
                raise ValueError("Failed to find polynomial")
    return (NbTilesPerHyp, MARSSpaceShiftPerHyp)


def generateSetMembershipTest(S):
    if S.is_empty():
        return "0"  # always false
    # Universe testing (empty test if so)
    # That's a two-way test : using ISL's coalescing mechanism, and decomposing the original set into basic sets
    # If either turns out an universe, we're good
    bsets = S.get_basic_sets()
    univ = [bs.is_universe() for bs in bsets]
    Scoal = S.copy().coalesce()
    if True in univ or Scoal.plain_is_universe():
        return "1"  # always true
    DomParamMap = (
        S.copy()
        .move_dims(
            dim_type.param, S.dim(dim_type.param), dim_type.set, 0, S.dim(dim_type.set)
        )
        .identity()
        .set_tuple_name(dim_type.in_, "T")
    )
    astBuild = AstBuild.alloc(DomParamMap.get_ctx())
    memberhsipTestStr = astBuild.ast_from_schedule(
        DomParamMap.to_union_map()
    ).to_C_str()
    m = re.match("if \((.*)\)", memberhsipTestStr)
    assert not m is None, "Test empty for non-universe or empty set"
    return m.group(1)


def getFlowOutMARSTransactionStart(FlowOutMARS):
    # Browse through all MARSes and each time a consumer is added, it marks a
    # new transaction start.
    PreviousConsumers = []
    NewTransactionIndices = []
    for index in range(len(FlowOutMARS)):
        (I, _) = FlowOutMARS[index]
        for C in I:  # C = consumer tile represented as list of crossed hyps
            if not C in PreviousConsumers:
                NewTransactionIndices += [index]
        PreviousConsumers = I
    return NewTransactionIndices


def generateUnrollDirectiveUnionMap(Schedule):
    # Schedule is wxpected to be a union map
    map_list = Schedule.get_map_list()
    unrollDirective = Map.from_range(Set("{ unroll[x] }"))
    codegenUMap = UnionMap("{ }")
    for i in range(map_list.n_map()):
        smap = map_list.get_map(i)
        smap_range_id = Map.from_domain(
            smap.range().get_space().universe_set()
        ).apply_range(unrollDirective)
        codegenUMap = codegenUMap.union(smap_range_id)
    return codegenUMap
