from mars.compute.utils import complement, scalarProduct, nonTrivialPartsOfSet
from mars.polyhedral.utils import getTileInnerSet

from islpy import Map, Set, dim_type

def getFlowOut(Space, Dependences, Hyperplanes, NormalVectors, TileSizes):
    # Flow-out per hyperplane
    FlowOuts = []
    SpaceStr = "[" + ", ".join(Space) + "]"
    for i in range(len(Hyperplanes)):
        H = Hyperplanes[i]
        c = NormalVectors[i]
        d = TileSizes[i]
        FlowOutHyp = []
        for b in Dependences:
            m = -1 * scalarProduct(c, b)
            if m != 0:
                while m < 0:
                    m += d

                Fij = Set(
                    "{ "
                    + SpaceStr
                    + " : ("
                    + H
                    + ") mod "
                    + str(d)
                    + " >= "
                    + str(m)
                    + " }"
                )
                FlowOutHyp.append(Fij)
            else:
                FlowOutHyp.append(Set("{ " + SpaceStr + " : false }"))
        FlowOuts.append(FlowOutHyp)

    TotalFlowOut = Set("{ " + SpaceStr + " : false }")
    for F in FlowOuts:
        for X in F:
            TotalFlowOut = TotalFlowOut.union(X.copy())

    return TotalFlowOut


def constrainDependenceToConsumerTile(
    b, T, Space, NormalVectors, Hyperplanes, TileSizes
):
    SpaceStr = "[" + ", ".join(Space) + "]"
    allHyperplanesCrossed = Set("{ " + SpaceStr + " }")
    for h in T:
        c = NormalVectors[h]
        H = Hyperplanes[h]
        d = TileSizes[h]
        m = -1 * scalarProduct(c, b)

        if m != 0:
            while m < 0:
                m += d
            allHyperplanesCrossed = allHyperplanesCrossed.intersect(
                Set(
                    "{ "
                    + SpaceStr
                    + " : (("
                    + H
                    + ") mod "
                    + str(d)
                    + " >= "
                    + str(m)
                    + ") }"
                )
            )
        else:
            allHyperplanesCrossed = allHyperplanesCrossed.intersect(
                Set("{ " + SpaceStr + " : false }")
            )
            return allHyperplanesCrossed  # There's nothing else to do here, let's go faster

    for h in range(len(Hyperplanes)):
        if h in T:
            continue
        c = NormalVectors[h]
        H = Hyperplanes[h]
        d = TileSizes[h]
        m = -1 * scalarProduct(c, b)

        if m != 0:
            while m < 0:
                m += d
            allHyperplanesCrossed = allHyperplanesCrossed.intersect(
                Set(
                    "{ "
                    + SpaceStr
                    + " : (not (("
                    + H
                    + ") mod "
                    + str(d)
                    + " >= "
                    + str(m)
                    + ")) }"
                )
            )

    return allHyperplanesCrossed


# Compute consumer tiles
def getConsumerTiles(Space, Dependences, Hyperplanes, NormalVectors, TileSizes):
    SpaceStr = "[" + ", ".join(Space) + "]"
    NeighboringTiles = nonTrivialPartsOfSet([i for i in range(len(Hyperplanes))])
    Consumers = []
    for T in NeighboringTiles:
        existsOneDependenceToT = Set("{ " + SpaceStr + " : false }")
        for b in Dependences:
            existsOneDependenceToT = existsOneDependenceToT.union(
                constrainDependenceToConsumerTile(
                    b, T, Space, NormalVectors, Hyperplanes, TileSizes
                )
            ).coalesce()
        if not existsOneDependenceToT.is_empty():
            # print("Consumer tile: " + str(T))
            Consumers.append(T)
    return Consumers


# Compute MARS, improved
def getMARS(Space, Dependences, Hyperplanes, NormalVectors, TileSizes):
    SpaceStr = "[" + ", ".join(Space) + "]"
    AllNeighboringTiles = nonTrivialPartsOfSet([i for i in range(len(Hyperplanes))])
    Consumers = getConsumerTiles(
        Space, Dependences, Hyperplanes, NormalVectors, TileSizes
    )

    def makeMARS(I, A):
        # I is the consumer tiles
        E = complement(I, AllNeighboringTiles)
        # ... and E is the non-consumer tiles

        existsNoDependenceToE = Set("{ " + SpaceStr + " }")
        for b in Dependences:
            bNoLeadToT = Set("{ " + SpaceStr + " }")
            for T in E:
                bNoLeadToT = bNoLeadToT.intersect(
                    constrainDependenceToConsumerTile(
                        b, T, Space, NormalVectors, Hyperplanes, TileSizes
                    ).complement()
                ).coalesce()
            existsNoDependenceToE = existsNoDependenceToE.intersect(
                bNoLeadToT
            ).coalesce()
        B = existsNoDependenceToE.copy()

        Res = A.intersect(B)

        return (I, Res)

    def getPointsConsumedAtLeastBy(I):
        A = Set("{ " + SpaceStr + " }")
        # For each consumer tile...
        for T in I:
            existsOneDependenceToT = Set("{ " + SpaceStr + " : false }")
            for b in Dependences:
                existsOneDependenceToT = existsOneDependenceToT.union(
                    constrainDependenceToConsumerTile(
                        b, T, Space, NormalVectors, Hyperplanes, TileSizes
                    )
                ).coalesce()
            A = A.intersect(existsOneDependenceToT).coalesce()

        return A

    def recurse(consumers, consumersFlowIn, others):
        allMARS = []
        thisConsumers = [T for T in consumers]
        queueForOthers = [T for T in others]
        # For each possible consumer tile not yet added to the list...
        while not len(queueForOthers) == 0:
            # Pop the added consumer and drop it from the list
            # Add the consumer tile to the list
            nCT = queueForOthers.pop(0)
            thisConsumers.insert(0, nCT)
            # print(thisConsumers)

            # Test if this set of consumer tiles has a common flow-in
            # commonFlowIn = getPointsConsumedAtLeastBy(thisConsumers)
            commonFlowIn = (
                consumersFlowIn.copy()
                .intersect(getPointsConsumedAtLeastBy([nCT]))
                .coalesce()
            )

            if not commonFlowIn.is_empty():
                # Make a MARS with exactly this set of consumers
                # (make an actual copy of the list)
                thisConsumersExactlyMARS = makeMARS(
                    [x for x in thisConsumers], commonFlowIn.copy()
                )
                if not thisConsumersExactlyMARS[1].is_empty():
                    allMARS += [thisConsumersExactlyMARS]

                # Carry on adding other consumers from those not yet added
                thisConsumersPlusOthersMARS = recurse(
                    thisConsumers, commonFlowIn.copy(), queueForOthers
                )
                allMARS += thisConsumersPlusOthersMARS
                # print("AllMars = " + str(allMARS))

            thisConsumers.pop(0)

        return allMARS

    return recurse([], Set("{ " + SpaceStr + " }"), Consumers)


# compute flow-in MARS
def getFlowInMARS(Space, Hyperplanes, ParamNames, TileSizes, MARS, ParamRelation=None):
    FlowInMARS = []
    for I, m in MARS:
        # I are the consumer tiles (list of crossed hyperplanes)
        # m is the MARS
        for T in I:
            # Go back to producer tile
            # GoBack = 0
            Decal = [0 for _ in Hyperplanes]
            # This assumes dependences cross hyperplanes forwards
            for i in T:
                Decal[i] = -1
            ProducerTileSpace = getTileInnerSet(
                Space, Hyperplanes, ParamNames, TileSizes, Decal
            )
            ThisMARSforThisProducer = m.copy().intersect(ProducerTileSpace).coalesce()
            if ParamRelation is None:
                ThisMARSforThisProducerGist = ThisMARSforThisProducer
            else:
                ThisMARSforThisProducerGist = ThisMARSforThisProducer.gist_params(
                    ParamRelation.copy()
                ).coalesce()
            if not ThisMARSforThisProducerGist.is_empty():
                FlowInMARS.append((I, ThisMARSforThisProducerGist, T))
    return FlowInMARS


def getFlowOutMARS(Space, Hyperplanes, ParamNames, TileSizes, MARS, ParamRelation=None):
    Tile = getTileInnerSet(
        Space, Hyperplanes, ParamNames, TileSizes, [0 for _ in Hyperplanes]
    )
    FlowOutMARS = []
    for I, m in MARS:
        ThisMARSforThisTile = m.copy().intersect(Tile.copy()).coalesce()
        if ParamRelation is None:
            ThisMARSforThisTileGist = ThisMARSforThisTile
        else:
            ThisMARSforThisTileGist = ThisMARSforThisTile.gist_params(
                ParamRelation.copy()
            ).coalesce()
        if not ThisMARSforThisTileGist.is_empty():
            FlowOutMARS.append((I, ThisMARSforThisTileGist))
    return FlowOutMARS


def annotateFlowInMARSWithRelations(FlowInMARSRel, Hyperplanes, ParamNames):
    NbParams = len(Hyperplanes)
    RelAsSet = {}
    for relation in FlowInMARSRel:
        RelAsSet[relation] = relation.copy().move_dims(
            dim_type.set, 0, dim_type.param, 0, NbParams
        )

    FlowInMARSRelAugmented = {}
    for relation in FlowInMARSRel:
        FlowInMARSAugmented = []
        FlowInMARS = FlowInMARSRel[relation]
        for I, M, T in FlowInMARS:
            # Compute the relation that the producer tile verifies
            Decal = [0 if i not in T else -1 for i in range(len(Hyperplanes))]
            DecalMap = Map(
                "{ ["
                + ", ".join(ParamNames)
                + "] -> ["
                + ", ".join(
                    [x + " - 1" if d == -1 else x for (x, d) in zip(ParamNames, Decal)]
                )
                + "] }"
            )
            RelComp = RelAsSet[relation].copy().apply(DecalMap)
            for R in FlowInMARSRel:
                if not RelComp.copy().intersect(RelAsSet[R].copy()).is_empty():
                    MatchingRel = R
                    break
            FlowInMARSAugmented += [(I, M, T, MatchingRel)]
        FlowInMARSRelAugmented[relation] = FlowInMARSAugmented
    return FlowInMARSRelAugmented
