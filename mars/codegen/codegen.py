import re
from islpy import Set, UnionMap, AstBuild, dim_type
from mars.codegen.utils import wrapIntoFunction, generateUnrollDirectiveUnionMap
from mars.polyhedral.barvinok import getBarvinokCardAsInt, getBarvinokCardAsIntWithGist


def generateMARSDataStructure(Dtype, FlowOutMARS):
    marsCard = [getBarvinokCardAsInt(M) for (_, M) in FlowOutMARS]
    marsDataStruct = []
    for m in range(len(FlowOutMARS)):
        c = marsCard[m]
        marsDataStruct += [f"   {Dtype} m{m}[{c}];"]
    return "struct mars {\n" + "\n".join(marsDataStruct) + "\n};"


def generateMARSDataStructureWithRelations(Dtype, FlowOutMARSRel, Relations):
    ret = ""
    for R in FlowOutMARSRel:
        rel_idx = Relations.index(R)
        FlowOutMARS = FlowOutMARSRel[R]
        marsCard = [getBarvinokCardAsIntWithGist(M, R) for (_, M) in FlowOutMARS]
        marsDataStruct = []
        for m in range(len(FlowOutMARS)):
            c = marsCard[m]
            marsDataStruct += [f"   {Dtype} m{m}[{c}];"]
        ret += (
            f"struct mars_R{rel_idx}"
            + " {\n"
            + "\n".join(marsDataStruct)
            + "\n} __attribute__((packed));\n"
        )
    return ret


def generateMARSGlobalAccessMacro(MARSSymbol, Hyperplanes, FlowOutMARS):
    nbParamsM = [M.dim(dim_type.param) for (_, M) in FlowOutMARS]
    assert min(nbParamsM) == max(nbParamsM), "Params differ by MARS"
    nbParams = max(nbParamsM)
    assert nbParams == len(Hyperplanes), "One parameter per tiling hyperplane"

    # Generate global MARS array access macro
    nbParams = max(nbParamsM)
    paramNames = [
        FlowOutMARS[0][1].get_dim_name(dim_type.param, i) for i in range(nbParams)
    ]
    spaceSize = [f"M{i+1}" for i in range(nbParams)]
    spaceLinearization = [
        "("
        + " * ".join(
            [spaceSize[j] for j in range(i + 1, nbParams)]
            if i < nbParams - 1
            else ["1"]
        )
        + " * ("
        + paramNames[i]
        + "))"
        for i in range(nbParams)
    ]

    spaceSizeMacros = "\n".join(
        [
            f"#define {S} 0 // Replace by number of tiles along hyperplane {H}"
            for (S, H) in list(zip(spaceSize, Hyperplanes))
        ]
    )
    marsMacro = (
        "#define marsGlobal("
        + ", ".join(paramNames)
        + ") ("
        + MARSSymbol
        + " + "
        + " + ".join(spaceLinearization)
        + ")"
    )
    return spaceSizeMacros + "\n" + marsMacro


def generateMARSRetrieveCode(
    FlowInMARS, FlowOutMARS, SourcePtr, DestArray, CounterVariable
):
    def tuplizeI(I):
        return tuple([tuple(x) for x in I])

    # Make transactions.

    FlowOutTuples = [tuplizeI(I) for (I, _) in FlowOutMARS]

    nbParams = FlowOutMARS[0][1].dim(dim_type.param)
    paramNames = [
        FlowOutMARS[0][1].get_dim_name(dim_type.param, i) for i in range(nbParams)
    ]

    producerTileChange = True
    isContiguous = True
    readLength = 0
    marsReadCode = ""
    subjectMARSNumbers = []
    for i in range(len(FlowInMARS)):
        (I, M, T) = FlowInMARS[i]
        Tdec = [0 for _ in range(nbParams)]
        for x in T:
            Tdec[x] = 1
        Tindices = [
            (paramNames[i] + "-" + str(Tdec[i])) if Tdec[i] != 0 else paramNames[i]
            for i in range(len(Tdec))
        ]

        if producerTileChange:
            marsReadCode += (
                f"{SourcePtr} = marsGlobal("
                + ", ".join(Tindices)
                + ")->m"
                + str(FlowOutTuples.index(tuplizeI(I)))
                + ";\n"
            )

        subjectMARSNumbers += [str(FlowOutTuples.index(tuplizeI(I)))]

        if i <= len(FlowInMARS) - 2:
            (Iy, _, Ty) = FlowInMARS[i + 1]
            if (
                T == Ty
                and FlowOutTuples.index(tuplizeI(I))
                == FlowOutTuples.index(tuplizeI(Iy)) - 1
            ):
                isContiguous = True
            else:
                isContiguous = False

            if not T == Ty:
                producerTileChange = True
            else:
                producerTileChange = False
        else:
            isContiguous = False

        readLength += getBarvinokCardAsInt(M)
        if not isContiguous:
            marsReadCode += (
                "/* Covered MARS: "
                + ", ".join(subjectMARSNumbers)
                + " */\n"
                + f"for(unsigned i = 0; i < {readLength}; ++i)"
                + "{\n"
                + f"  {DestArray}[{CounterVariable}] = *{SourcePtr};\n"
                + f"  {CounterVariable}++;\n"
                + f"  {SourcePtr}++;\n"
                + "}\n"
            )
            readLength = 0
            subjectMARSNumbers = []

    return marsReadCode


def generateMARSStoreCode(
    FlowOutMARS, DestPtr, SourceArray, CounterVariable, ConvertMARSGlobalType=None
):
    nbParams = FlowOutMARS[0][1].dim(dim_type.param)
    paramNames = [
        FlowOutMARS[0][1].get_dim_name(dim_type.param, i) for i in range(nbParams)
    ]

    Conversion = ""
    if ConvertMARSGlobalType is not None:
        Conversion = f"({ConvertMARSGlobalType})"

    totalCard = sum([getBarvinokCardAsInt(M) for (_, M) in FlowOutMARS])
    marsWriteCode = (
        f"{DestPtr} = {Conversion}(marsGlobal("
        + ", ".join(paramNames)
        + "));\n"
        + f"for(unsigned i = 0; i < {totalCard}; ++i)"
        + "{\n"
        + f"  *{DestPtr} = {SourceArray}[{CounterVariable}];\n"
        + f"  {CounterVariable}++;\n"
        + f"  {DestPtr}++;\n"
        + "}\n"
    )

    return marsWriteCode


def generateMARSInputCode(
    Dtype,
    StatementWriteUnionMap,
    MARSSymbol,
    FlowInMARS,
    FlowOutMARS,
    Unroll=False,
    SuplRel=[],
    LocalArrayBounds={},
    TestForInset=False,
    UnrollIntoROM=False,
    UseGlobalPointer=True,
):
    ParamSpace = Set.universe(FlowInMARS[0][1].params().get_space())
    Relations = [ParamSpace]
    FlowInMARSRel = {ParamSpace: [(M, T, P, ParamSpace) for (M, T, P) in FlowInMARS]}
    FlowOutMARSRel = {ParamSpace: FlowOutMARS}
    return generateMARSInputCodeWithRelations(
        Dtype,
        StatementWriteUnionMap,
        MARSSymbol,
        FlowInMARSRel,
        FlowOutMARSRel,
        Relations,
        {},
        ParamSpace,
        Unroll=Unroll,
        SuplRelForeach=SuplRel,
        LocalArrayBounds=LocalArrayBounds,
        TestForInset=TestForInset,
        UnrollIntoROM=UnrollIntoROM,
        UseGlobalPointer=UseGlobalPointer,
    )


def generateMARSInputFunction(
    Dtype,
    StatementWriteUnionMap,
    MARSSymbol,
    FlowInMARS,
    FlowOutMARS,
    unroll=False,
    supl_rel=[],
    local_array_size_min_max={},
    test_for_inset=False,
):
    return wrapIntoFunction(
        generateMARSInputCode(
            Dtype,
            StatementWriteUnionMap,
            MARSSymbol,
            FlowInMARS,
            FlowOutMARS,
            Unroll=unroll,
            SuplRel=supl_rel,
            getLocalArrayBounds=local_array_size_min_max,
            TestForInset=test_for_inset,
        )
    )


def generateMARSInputCodeWithRelations(
    Dtype,
    StatementWriteUnionMap,
    MARSSymbol,
    FlowInMARSRel,
    FlowOutMARSRel,
    Relations,
    RelationNames,
    SelectedRelation,
    Unroll=False,
    SuplRelForeach=[],
    LocalArrayBounds={},
    TestForInset=False,
    UnrollIntoROM=False,
    UseGlobalPointer=True,
):
    def tuplizeI(I):
        return tuple([tuple(x) for x in I])

    # Check compatibility of params
    assert not (
        UnrollIntoROM and not Unroll
    ), "Unroll must be enabled to use ROM MARS address storage"

    # Check that all MARS have the same number of params
    nbParamsM = [
        M.dim(dim_type.param) for R in FlowOutMARSRel for (_, M) in FlowOutMARSRel[R]
    ]
    assert min(nbParamsM) == max(nbParamsM), "Params differ by MARS"

    # Get parameter names
    nbParams = max(nbParamsM)
    paramNames = [
        FlowOutMARSRel[Relations[0]][0][1].get_dim_name(dim_type.param, i)
        for i in range(nbParams)
    ]
    functionArgsArr = ["mars", "*" + MARSSymbol] + [("int", P) for P in paramNames]

    # Get symbols of on-chip memories
    innerSymbols = []
    symbolDims = []
    writeMaps = StatementWriteUnionMap.get_map_list()
    for i in range(writeMaps.n_map()):
        writeMap = writeMaps.get_map(i)
        TargetSymbol = writeMap.get_tuple_name(dim_type.out)
        if not TargetSymbol in innerSymbols:
            innerSymbols += [TargetSymbol]
            symbolDims += [writeMap.dim(dim_type.out)]

    # Make a constraint map to order them in alphabetical order
    sortedSymbols = [x for x in innerSymbols]
    sortedSymbols.sort()
    symbolAlphaConstraints = []
    for i in range(len(sortedSymbols)):
        for j in range(i + 1, len(sortedSymbols)):
            dimsI = ", ".join([f"i{x}" for x in range(symbolDims[i])])
            dimsJ = ", ".join(
                [f"i{x}" for x in range(symbolDims[i], symbolDims[i] + symbolDims[j])]
            )
            symbolAlphaConstraints += [
                f"{sortedSymbols[i]}[{dimsI}] -> {sortedSymbols[j]}[{dimsJ}]"
            ]

    symbolAlphaConstraintMap = UnionMap("{ " + "; ".join(symbolAlphaConstraints) + " }")

    # Build function arguments
    functionArgsArr += [
        (
            Dtype,
            innerSymbols[i]
            + (
                "[]" * symbolDims[i]
                if not innerSymbols[i] in LocalArrayBounds
                else "".join(
                    ["[" + str(d) + "]" for d in LocalArrayBounds[innerSymbols[i]][0]]
                )
            ),
        )
        for i in range(len(innerSymbols))
    ]

    # Make tuplized lists
    FlowOutTuples = {
        R: [tuplizeI(I) for (I, _) in FlowOutMARSRel[R]] for R in FlowOutMARSRel
    }
    MARSCard = {
        R: [getBarvinokCardAsIntWithGist(M, R) for (_, M) in FlowOutMARSRel[R]]
        for R in FlowOutMARSRel
    }

    if len(SuplRelForeach) == 0:
        x = Set.universe(Relations[0].get_space())
        SuplRelForeach = [x]
    # for R in Relations:

    CopyStmts = []
    ROMFillIn = ""
    CurrentROMPos = 0
    # Get consumer relation name
    ConsRelationName = (
        ""
        if not SelectedRelation in RelationNames
        else f"_{RelationNames[SelectedRelation]}"
    )
    # Get flow-in MARS for this relation
    FlowInMARS = FlowInMARSRel[SelectedRelation]
    for Rs in SuplRelForeach:
        # Is there any supplemental restriction?
        rs_is_empty = (
            Rs.plain_is_universe()
        )  # Emptiness here means "emptiness of conditions"

        # Generate code to browse the MARS
        Ctx = FlowInMARS[0][1].get_ctx()
        islAstBuild = AstBuild.alloc(
            Ctx
        )  # .set_options(UnionMap(codegenOptionsUnionMap))
        RsCopyStmts = []
        ROMFillInStmts = ""
        isContiguous = False

        if not rs_is_empty:
            # Add info as comment
            RsCopyStmts += ["    /* Supplemental Relation: " + str(Rs) + " */\n"]
            # Make guards
            Nparams_Rs = Rs.dim(dim_type.param)
            Rs_map = UnionMap.from_range(
                Rs.copy()
                .move_dims(dim_type.set, 0, dim_type.param, 0, Nparams_Rs)
                .identity()
                .move_dims(dim_type.param, 0, dim_type.in_, 0, Nparams_Rs)
                .range()
            )
            Guard = islAstBuild.ast_from_schedule(Rs_map)

        # Compute MARS read-in lengths
        CoalescedReadLength = []
        CurrentLength = 0
        for i in range(len(FlowInMARS)):
            if not isContiguous and CurrentLength > 0:
                CoalescedReadLength += [CurrentLength]
                CurrentLength = 0
            isContiguous = False
            (Ix, _, Tx, Rx) = FlowInMARS[i]
            CurrentLength += MARSCard[Rx][FlowOutTuples[Rx].index(tuplizeI(Ix))]
            if i <= len(FlowInMARS) - 2:
                (Iy, _, Ty, Ry) = FlowInMARS[i + 1]
                if (
                    Tx == Ty
                    and Rx == Ry
                    and FlowOutTuples[Rx].index(tuplizeI(Ix))
                    == FlowOutTuples[Rx].index(tuplizeI(Iy)) - 1
                ):
                    isContiguous = True
        # Add last transaction
        CoalescedReadLength += [CurrentLength]

        # Use flow-in MARS in provided order (e.g., from LP)
        CurrentTransaction = 0
        producerTileChange = True

        for i in range(len(FlowInMARS)):
            (I, M, T, Rprod) = FlowInMARS[i]

            thisMARSCode = ""
            # Reset pointer if not contiguous
            if not isContiguous:
                Tdec = [0 for i in range(len(paramNames))]
                for x in T:
                    Tdec[x] = 1
                Tindices = [
                    (paramNames[i] + "-" + str(Tdec[i]))
                    if Tdec[i] != 0
                    else paramNames[i]
                    for i in range(len(Tdec))
                ]

                # Shall we generate an inset guard?
                thisMARSCode = ""
                if TestForInset and producerTileChange:
                    thisMARSCode += (
                        "  if(PRODUCER_IS_INSET(" + ", ".join(Tindices) + ")) {\n"
                    )

                RelationName = (
                    "" if not Rprod in RelationNames else f"_{RelationNames[Rprod]}"
                )

                if not UnrollIntoROM:
                    if UseGlobalPointer:
                        thisMARSCode += (
                            f"  marsPtr = marsGlobal{RelationName}("
                            + ", ".join(Tindices)
                            + ")->m"
                            + str(FlowOutTuples[Rprod].index(tuplizeI(I)))
                            + ";"
                        )
                else:
                    ROMFillInStmts += "\n"
                    if UseGlobalPointer:
                        ROMFillInStmts += (
                            f"marsPtr = marsGlobal{RelationName}("
                            + ", ".join(Tindices)
                            + ")->m"
                            + str(FlowOutTuples[Rprod].index(tuplizeI(I)))
                            + ";\n"
                        )
                    ROMFillInStmts += (
                        f"for(int i={CurrentROMPos}; i<{CurrentROMPos+CoalescedReadLength[CurrentTransaction]}; ++i) {{\n"
                        + f"    struct mars_transfert mt = FPGA_MARS_IN_TBL{ConsRelationName}[i];\n"
                        + (
                            f"    {Dtype} val = *marsPtr; \n"
                            if UseGlobalPointer
                            else ""
                        )
                        + (
                            (
                                "    switch (mt.array) {\n"
                                + "".join(
                                    f"        case MARS_DATA_ENUM::{name}: {{\n"
                                    + f"            marsToMem_{name}("
                                    + ", ".join([f"mt.dim{i}" for i in range(dims)])
                                    + (", val);\n" if UseGlobalPointer else ");\n")
                                    + "            break;\n"
                                    + "        }\n"
                                    for name, dims in zip(innerSymbols, symbolDims)
                                )
                                + "    }\n"
                            )
                            if len(innerSymbols) > 1
                            else f"    marsToMem_{list(innerSymbols)[0]}("
                            + ", ".join([f"mt.dim{i}" for i in range(symbolDims[0])])
                            + (", val);\n" if UseGlobalPointer else ");\n")
                        )
                        + ("    marsPtr++;\n" if UseGlobalPointer else "")
                        + "}"
                    )
                    CurrentROMPos += CoalescedReadLength[CurrentTransaction]
                thisMARSCode += (
                    "\n/* This transaction will be "
                    + str(CoalescedReadLength[CurrentTransaction])
                    + " words long */\n"
                )
                CurrentTransaction += 1
            isContiguous = False

            # Build access code
            ParamRelComb = SelectedRelation.copy().intersect(Rs.copy())
            TargetCells = (
                M.copy()
                .apply(StatementWriteUnionMap.copy())
                .gist_params(ParamRelComb)
                .coalesce()
            )
            assert getBarvinokCardAsIntWithGist(
                TargetCells, ParamRelComb
            ) == getBarvinokCardAsIntWithGist(
                M, ParamRelComb
            ), "Illegal local allocation (overwriting)"

            # Order the points inside MARSes according to alphabetical and
            # lexicographic order

            Schedule = (
                TargetCells.compute_schedule(symbolAlphaConstraintMap, UnionMap("{ }"))
                .get_root()
                .get_subtree_schedule_union_map()
            )
            if Unroll:
                codegenOptionsUnionMap = generateUnrollDirectiveUnionMap(Schedule)
            else:
                codegenOptionsUnionMap = UnionMap("{ }")
            scheduleMaps = Schedule.get_map_list()
            newSchedule = UnionMap.empty_ctx(Schedule.get_ctx())
            for j in range(scheduleMaps.n_map()):
                thisMap = scheduleMaps.get_map(j)
                thisMapMod = thisMap
                if not UnrollIntoROM:
                    TargetArray = thisMap.get_tuple_name(dim_type.in_)
                    thisMapMod = thisMap.set_tuple_name(
                        dim_type.in_, "marsToMem_" + TargetArray
                    )
                newSchedule = newSchedule.add_map(thisMapMod)

            islAstBuilderForThisMARS = AstBuild.from_context(ParamRelComb).set_options(
                codegenOptionsUnionMap
            )

            # Build an AST describing a MARS
            # It calls marsToMem as a macro that will perform the actual copy
            AST = islAstBuilderForThisMARS.ast_from_schedule(newSchedule)
            if not UnrollIntoROM:
                # The AST is to be executed straight.
                thisMARSCode += "  " + AST.to_C_str().replace("\n", "\n  ")[:-2]
            else:
                # The AST is turned into a set of array cells.
                new_str = (
                    AST.to_C_str()
                    .replace("{\n", "")
                    .replace("}\n", "")
                    .replace("  \n", "")
                    .replace("  ", "")
                    .replace("\n", " ")
                    .replace("(", "{")
                    .replace(")", "}")
                    .replace(";", ",")
                )
                thisMARSCode += new_str
            # Check for contiguity
            if i <= len(FlowInMARS) - 2:
                (Iy, _, Ty, Ry) = FlowInMARS[i + 1]
                if (
                    T == Ty
                    and Rprod == Ry
                    and FlowOutTuples[Rprod].index(tuplizeI(I))
                    == FlowOutTuples[Rprod].index(tuplizeI(Iy)) - 1
                ):
                    isContiguous = True
                if not T == Ty:
                    producerTileChange = True
                    if TestForInset:
                        thisMARSCode += "  }\n"
                else:
                    producerTileChange = False
            elif TestForInset:
                thisMARSCode += "  }\n"

            RsCopyStmts += [thisMARSCode]

        if not rs_is_empty:
            if UnrollIntoROM:
                CopyStmts += RsCopyStmts
                ROMFillIn += "  " + Guard.to_C_str().replace(
                    "  ();",
                    "  {\n" + ROMFillInStmts.replace("\n", "\n  ") + "\n  }\n",
                )
            else:
                CopyStmts += [
                    "  "
                    + Guard.to_C_str().replace(
                        "  ();",
                        "  {\n"
                        + "".join(RsCopyStmts).replace("\n", "\n  ")
                        + "\n  }\n",
                    )
                ]
        else:
            CopyStmts += RsCopyStmts
            if UnrollIntoROM:
                ROMFillIn += ROMFillInStmts

    MARSPointerDecl = (Dtype + "* marsPtr;\n") if UseGlobalPointer else ""
    if not UnrollIntoROM:
        return (MARSPointerDecl + "".join(CopyStmts), functionArgsArr)
    else:
        return (
            re.sub(
                "([a-zA-Z0-9-_]+){",
                "{MARS_DATA_ENUM::\\1, ",
                "".join(CopyStmts)[:-2],
            ),
            functionArgsArr,
            MARSPointerDecl + ROMFillIn,
        )


def generateMARSOutputCode(
    Dtype,
    StatementWriteUnionMap,
    MARSSymbol,
    FlowOutMARS,
    Unroll=False,
    SuplRel=[],
    LocalArrayBounds={},
    MarsGlobalSuffix="",
    AllRel=None,
    UnrollIntoROM=False,
    UseGlobalPointer=True,
):
    def tuplizeI(I):
        return tuple([tuple(x) for x in I])

    # Make tuplized lists
    # FlowOutTuples = [tuplizeI(I) for (I, _) in FlowOutMARS]
    if AllRel is None:
        MARSCard = [getBarvinokCardAsInt(M) for (_, M) in FlowOutMARS]
    else:
        MARSCard = [
            getBarvinokCardAsInt(M.gist_params(AllRel)) for (_, M) in FlowOutMARS
        ]

    # Check that all MARS have the same number of params
    nbParamsM = [M.dim(dim_type.param) for (_, M) in FlowOutMARS]
    assert min(nbParamsM) == max(nbParamsM), "Params differ by MARS"

    # Get parameter names
    nbParams = max(nbParamsM)
    paramNames = [
        FlowOutMARS[0][1].get_dim_name(dim_type.param, i) for i in range(nbParams)
    ]
    functionArgsArr = [("struct mars", "*" + MARSSymbol)] + [
        ("int", P) for P in paramNames
    ]

    # Get used (on-chip) symbols and their dimensions
    innerSymbols = []
    symbolDims = []
    writeMaps = StatementWriteUnionMap.get_map_list()
    allTargetSymbols = set()
    for i in range(writeMaps.n_map()):
        writeMap = writeMaps.get_map(i)
        TargetSymbol = writeMap.get_tuple_name(dim_type.out)
        allTargetSymbols.add(TargetSymbol)
        if not TargetSymbol in innerSymbols:
            innerSymbols += [TargetSymbol]
            symbolDims += [writeMap.dim(dim_type.out)]

    # Make a constraint map to order them in alphabetical order
    sortedSymbols = [x for x in innerSymbols]
    sortedSymbols.sort()
    symbolAlphaConstraints = []
    for i in range(len(sortedSymbols)):
        for j in range(i + 1, len(sortedSymbols)):
            dimsI = ", ".join([f"i{x}" for x in range(symbolDims[i])])
            dimsJ = ", ".join(
                [f"i{x}" for x in range(symbolDims[i], symbolDims[i] + symbolDims[j])]
            )
            symbolAlphaConstraints += [
                f"{sortedSymbols[i]}[{dimsI}] -> {sortedSymbols[j]}[{dimsJ}]"
            ]

    symbolAlphaConstraintMap = UnionMap("{ " + "; ".join(symbolAlphaConstraints) + " }")

    # Function arguments
    functionArgsArr += [
        (
            Dtype,
            innerSymbols[i]
            + (
                "[]" * symbolDims[i]
                if not innerSymbols[i] in LocalArrayBounds
                else "".join(
                    ["[" + str(d) + "]" for d in LocalArrayBounds[innerSymbols[i]][0]]
                )
            ),
        )
        for i in range(len(innerSymbols))
    ]

    # If there is no given relation between the parameters, then let's just add
    # a virtual one (universe means... no constraint at all on parameters)
    if len(SuplRel) == 0:
        SuplRel = [Set.universe(FlowOutMARS[0][1].params().get_space())]

    # No relation means we'll intersect parameter spaces with the N-dimensional universe
    # That way we'll not alter MARS if there is no relation between params to be verified
    if AllRel is None:
        AllRel = Set.universe(FlowOutMARS[0][1].params().get_space())

    CopyStmts = []
    # Intialize MARS pointer
    if MarsGlobalSuffix != "":
        MarsGlobalSuffix = "_" + MarsGlobalSuffix
    if not UnrollIntoROM and UseGlobalPointer:
        CopyStmts += [
            "  marsPtr = ("
            + Dtype
            + "*)(marsGlobal"
            + MarsGlobalSuffix
            + "("
            + ", ".join(paramNames)
            + ")); /* This transaction will be "
            + str(sum(MARSCard))
            + " words long */\n"
        ]
    else:
        CopyStmts += [
            "\n/* This transaction will be " + str(sum(MARSCard)) + " words long */\n"
        ]

    ROMFillIn = ""
    CurrentROMPos = 0
    MARSTotalCard = sum(MARSCard)
    # Browse MARS with each of the given relations
    for Rs in SuplRel:
        # Is there any supplemental restriction?
        rs_is_empty = (
            Rs.plain_is_universe()
        )  # Emptiness here means "emptiness of conditions"
        RsCopyStmts = []

        # Generate code to browse the MARS
        Ctx = FlowOutMARS[0][1].get_ctx()
        islAstBuild = AstBuild.alloc(
            Ctx
        )  # .set_options(UnionMap(codegenOptionsUnionMap))

        if not rs_is_empty:
            # Add info as comment
            RsCopyStmts += ["    /* Supplemental Relation: " + str(Rs) + " */\n"]
            # Make guards
            Nparams_Rs = Rs.dim(dim_type.param)
            Rs_map = UnionMap.from_range(
                Rs.copy()
                .move_dims(dim_type.set, 0, dim_type.param, 0, Nparams_Rs)
                .identity()
                .move_dims(dim_type.param, 0, dim_type.in_, 0, Nparams_Rs)
                .range()
            )
            Guard = islAstBuild.ast_from_schedule(Rs_map)

        # Generate ROM Fill-in code for this relation
        ROMFillInStmts = (
            "\n"
            + (
                Dtype
                + " *marsPtr = ("
                + Dtype
                + "*)(marsGlobal"
                + MarsGlobalSuffix
                + "("
                + ", ".join(paramNames)
                + "));\n"
                if UseGlobalPointer
                else ""
            )
            + (
                f"for(int i={CurrentROMPos}; i<{CurrentROMPos+MARSTotalCard}; ++i) {{\n"
                + f"    struct mars_transfert mt = FPGA_MARS_OUT_TBL{MarsGlobalSuffix}[i];\n"
                + (f"    {Dtype} val; \n" if UseGlobalPointer else "")
                + (
                    (
                        "    switch (mt.array) {\n"
                        + "".join(
                            f"        case MARS_DATA_ENUM::{name}: {{\n"
                            + f"            memToMARS_{name}("
                            + ", ".join([f"mt.dim{i}" for i in range(dims)])
                            + (", val);\n" if UseGlobalPointer else ");\n")
                            + "            break;\n"
                            + "        }\n"
                            for name, dims in zip(innerSymbols, symbolDims)
                        )
                        + "    }\n"
                    )
                    if len(allTargetSymbols) > 1
                    else f"    memToMARS_{innerSymbols[0]}("
                    + ", ".join([f"mt.dim{i}" for i in range(symbolDims[0])])
                    + (", val);\n" if UseGlobalPointer else ");\n")
                )
                + ("    *marsPtr = val;\n" if UseGlobalPointer else "")
                + ("    marsPtr++;\n" if UseGlobalPointer else "")
                + "}"
            )
        )
        CurrentROMPos += MARSTotalCard

        # Use flow-out MARS in provided order (from LP)
        for _, M in FlowOutMARS:
            # Build access code
            TargetCells = (
                M.copy()
                .apply(StatementWriteUnionMap.copy())
                .gist_params(Rs.copy().intersect(AllRel.copy()))
                .coalesce()
            )
            assert getBarvinokCardAsIntWithGist(
                TargetCells, Rs.intersect(AllRel.copy())
            ) == getBarvinokCardAsIntWithGist(
                M, Rs.intersect(AllRel.copy())
            ), "Illegal local allocation (redundancy)"
            Schedule = (
                TargetCells.compute_schedule(symbolAlphaConstraintMap, UnionMap("{ }"))
                .get_root()
                .get_subtree_schedule_union_map()
            )
            if Unroll:
                codegenOptionsUnionMap = generateUnrollDirectiveUnionMap(Schedule)
            else:
                codegenOptionsUnionMap = UnionMap("{ }")
            scheduleMaps = Schedule.get_map_list()
            newSchedule = UnionMap.empty_ctx(Schedule.get_ctx())
            for j in range(scheduleMaps.n_map()):
                thisMap = scheduleMaps.get_map(j)
                thisMapMod = thisMap
                if not UnrollIntoROM:
                    TargetArray = thisMap.get_tuple_name(dim_type.in_)
                    thisMapMod = thisMap.set_tuple_name(
                        dim_type.in_, "memToMARS_" + TargetArray
                    )
                newSchedule = newSchedule.add_map(thisMapMod)

            islAstBuilderForThisMARS = AstBuild.from_context(Rs).set_options(
                codegenOptionsUnionMap
            )
            AST = islAstBuilderForThisMARS.ast_from_schedule(newSchedule)
            if not UnrollIntoROM:
                thisMARSCode = "  " + AST.to_C_str().replace("\n", "\n  ")[:-2]
            else:
                new_str = (
                    AST.to_C_str()
                    .replace("{\n", "")
                    .replace("}\n", "")
                    .replace("  \n", "")
                    .replace("  ", "")
                    .replace("\n", " ")
                    .replace("(", "{")
                    .replace(")", "}")
                    .replace(";", ",")
                )
                thisMARSCode = new_str
            RsCopyStmts += [thisMARSCode]

        # If there is a relation, then generate special code for it
        if not rs_is_empty:
            if UnrollIntoROM:
                CopyStmts += RsCopyStmts
                ROMFillIn += "  " + Guard.to_C_str().replace(
                    "  ();",
                    "  {\n" + ROMFillInStmts.replace("\n", "\n  ") + "\n  }\n",
                )
            else:
                CopyStmts += [
                    "  "
                    + Guard.to_C_str().replace(
                        "  ();",
                        "  {\n"
                        + "".join(RsCopyStmts).replace("\n", "\n  ")
                        + "\n  }\n",
                    )
                ]
        else:
            CopyStmts += RsCopyStmts
            if UnrollIntoROM:
                ROMFillIn += ROMFillInStmts

    if not UnrollIntoROM:
        return (Dtype + "* marsPtr;\n" + "".join(CopyStmts), functionArgsArr)
    else:
        return (
            re.sub(
                "([a-zA-Z0-9-_]+){",
                "{MARS_DATA_ENUM::\\1, ",
                "".join(CopyStmts)[:-2],
            ),
            functionArgsArr,
            ROMFillIn,
        )


def generateMARSOutputFunction(
    Dtype,
    StatementWriteUnionMap,
    MARSSymbol,
    FlowOutMARS,
    Unroll=False,
    SuplRel=[],
    getLocalArrayBounds={},
    CheckBounds=False,
):
    return wrapIntoFunction(
        generateMARSOutputCode(
            Dtype,
            StatementWriteUnionMap,
            MARSSymbol,
            FlowOutMARS,
            Unroll=Unroll,
            SuplRel=SuplRel,
            LocalArrayBounds=getLocalArrayBounds,
            CheckBounds=CheckBounds,
        )
    )
