from islpy import Map, Set, UnionMap, UnionSet, dim_type
from simplejson.tests import test_tuple
from mars.compute.utils import scalarProduct

def checkTilingLegal(Dependences, NormalVectors):
    AllTilingLegal = True
    for n in NormalVectors:
        allPos = True
        allNeg = False
        for b in Dependences:
            m = scalarProduct(n, b)
            if m > 0:
                allNeg = False
            if m < 0:
                allPos = False
        OK = allNeg or allPos
        if not OK:
            print("Tiling hyperplane with normal vector " + str(n) + " is illegal")
            AllTilingLegal = False
    #        else:
    #            print("Tiling hyperplane with normal vector " + str(n) + " is legal")
    return AllTilingLegal

def checkMARSEquality(A, B):
    if len(A) != len(B):
        return False

    for IA, MA in A:
        foundEquivalent = False
        for IB, MB in B:
            if len(IB) != len(IA):
                continue

            # Check that consumer tiles are the same
            for T in IB:
                if not T in IA:
                    continue

            # Check that MARS spaces are equal
            if not MB.is_equal(MA):
                continue

            # If we reached there, we have found an equivalent
            foundEquivalent = True
            break

        if not foundEquivalent:
            return False
    return True


def getDependencesAsUnionMap(Space, Dependences):
    M = UnionMap("{  : false }")
    for b in Dependences:
        Expr = [
            idx + " + " + str(bx)
            if bx > 0
            else idx + " - " + str(-bx)
            if bx < 0
            else idx
            for (idx, bx) in zip(Space, b)
        ]
        M = M.add_map(
            Map("{ [" + ", ".join(Space) + "] -> [" + ", ".join(Expr) + "] }")
        )
    return M


def getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal):
    SpaceStr = "[" + ", ".join(Space) + "]"
    TileAsParam = []
    for i in range(len(Hyperplanes)):
        H = Hyperplanes[i]
        d = TileSizes[i]
        TileAsParam.append(
            str(d)
            + " * ("
            + ParamNames[i]
            + " + "
            + str(Decal[i])
            + ") <= ("
            + H
            + ") < "
            + str(d)
            + " * (1 + "
            + ParamNames[i]
            + " + "
            + str(Decal[i])
            + ")"
        )

    ParamTile = (
        "["
        + ", ".join(ParamNames)
        + "] -> { "
        + SpaceStr
        + " : "
        + " and ".join(TileAsParam)
        + " }"
    )
    return Set(ParamTile)


def getLocalArrayBounds(
    Space,
    Hyperplanes,
    ParamNames,
    TileSizes,
    ReadMap,
    WriteMap,
    StmtRemap=None,
    gist_rel=None,
):
    Decal = [0 for _ in range(len(Hyperplanes))]
    # Get a tile in the original space, translated if need be
    Tile = getTileInnerSet(Space, Hyperplanes, ParamNames, TileSizes, Decal)
    if not gist_rel is None:
        Tile = Tile.gist_params(gist_rel)
    if not StmtRemap is None:
        Tile = Tile.apply(StmtRemap)

    # Compute the footprint of the tile
    Footprint = (
        Tile.copy().apply(WriteMap).union(Tile.copy().apply(ReadMap))
    )  # .drop_unused_params() cannot be applied to union sets

    symbol_bounds = {}

    # For each symbol...
    for i in range(Footprint.n_set()):
        s = Footprint.get_set_list().get_set(i)
        symbol = s.get_tuple_name()
        Ndim = s.dim(dim_type.set)
        dim_bounds = []
        # For each dimension, get min/max
        for j in range(s.dim(dim_type.set)):
            if j > 0:
                sproj = s.copy().project_out(dim_type.set, 0, j)
            else:
                sproj = s
            if j < Ndim - 1:
                sproj = sproj.project_out(dim_type.set, 1, Ndim - j - 1)
            lmin = sproj.lexmin().drop_unused_params()
            lmax = sproj.lexmax().drop_unused_params()
            min_aff_m = lmin.as_pw_multi_aff()
            max_aff_m = lmax.as_pw_multi_aff()
            min_aff = min_aff_m.as_multi_aff()
            max_aff = max_aff_m.as_multi_aff()

            if symbol in symbol_bounds:
                (cur_min_aff, cur_max_aff) = symbol_bounds[symbol][j]
                if min_aff.plain_cmp(cur_min_aff) < 0:
                    cur_min_aff = min_aff
                if max_aff.plain_cmp(cur_max_aff) > 0:
                    cur_max_aff = max_aff
                dim_bounds += [(cur_min_aff, cur_max_aff)]
            else:
                dim_bounds += [(min_aff, max_aff)]
        symbol_bounds[symbol] = dim_bounds

    symbol_sizeminmax = {}
    # Compute differences
    for symbol in symbol_bounds:
        DimExtrema = symbol_bounds[symbol]
        DimSizes = []
        DimMin = []
        DimMax = []
        for dmin, dmax in DimExtrema:
            assert not (dmin.involves_dims(dim_type.param, 0, dmin.dim(dim_type.param)))
            assert not (dmax.involves_dims(dim_type.param, 0, dmax.dim(dim_type.param)))
            dmin_val = dmin.get_constant_multi_val().get_val(0).to_python()
            dmax_val = dmax.get_constant_multi_val().get_val(0).to_python()
            diff = (
                1
                + dmax.as_pw_multi_aff()
                .sub(dmin.as_pw_multi_aff())
                .as_multi_aff()
                .get_constant_multi_val()
                .get_val(0)
                .to_python()
            )
            DimSizes += [diff]
            DimMin += [dmin_val]
            DimMax += [dmax_val]
        symbol_sizeminmax[symbol] = (DimSizes, DimMin, DimMax)
    return symbol_sizeminmax


# Compute tile origin as an affine map.
# We'll subtract this from every tile so that the FPGA tile loop nest
# is as simple as can be.
def getTileOriginTranslationFromHypInter(Space, Hyperplanes, TileSizes, ParamNames):
    SpaceStr = "[" + ", ".join(Space) + "]"
    TileAsParam = []
    for i in range(len(Hyperplanes)):
        H = Hyperplanes[i]
        d = TileSizes[i]
        TileAsParam.append(str(d) + " * (" + ParamNames[i] + ") = (" + H + ")")

    ParamTile = (
        "["
        + ", ".join(ParamNames)
        + "] -> { "
        + SpaceStr
        + " : "
        + " and ".join(TileAsParam)
        + " }"
    )
    Transl = Set(ParamTile).translation().reverse()
    return Transl


def getTileOriginTranslationFromLexmin(Space, Hyperplanes, TileSizes, ParamNames):
    oneTileLexmin = (
        getTileInnerSet(
            Space, Hyperplanes, ParamNames, TileSizes, [0 for _ in ParamNames]
        )
        .lexmin()
        .translation()
        .reverse()
    )
    return oneTileLexmin


def makeConstantTranslation(M, ConstantShift):
    ConstantTransl = (
        Set("{ [" + ", ".join([str(i) for i in ConstantShift]) + "] }")
        .translation()
        .reverse()
    )
    return M.apply_range(ConstantTransl)


def getSetFromUnionWithTupleName(S, Tuple):
    sets = [S.get_set_list().get_at(i) for i in range(S.n_set())]
    RUS = UnionSet.empty(S.get_space())
    for s in sets:
        if s.get_tuple_name() == Tuple:
            RUS = RUS.union(s)
    return Set.from_union_set(RUS)
