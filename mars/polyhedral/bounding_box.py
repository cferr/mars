from islpy import Set, UnionSet, dim_type


def getBoundingBox(Footprint, Symbol):
    fp = UnionSet.empty(Footprint.get_space())
    footprint_bmaps = [
        Footprint.get_set_list().get_at(i) for i in range(Footprint.n_set())
    ]
    for x in footprint_bmaps:
        if x.get_tuple_name() == Symbol:
            fp = fp.union(x)
    fpsym = Set.from_union_set(fp)
    bbsym = Set.universe(fpsym.get_space())
    ndimsym = fpsym.dim(dim_type.set)
    for i in range(ndimsym):
        fpdimmap = Set.universe(fpsym.get_space()).identity()
        # project out all dims except i
        if i > 0:
            fpsymdim = fpsym.project_out(dim_type.set, 0, i)
            fpdimmap = fpdimmap.project_out(dim_type.out, 0, i)
        else:
            fpsymdim = fpsym

        if i < ndimsym - 1:
            fpsymdim = fpsymdim.project_out(dim_type.set, 1, ndimsym - i - 1)
            fpdimmap = fpdimmap.project_out(dim_type.out, 1, ndimsym - i - 1)

        dimmin = fpsymdim.lexmin()
        dimmax = fpsymdim.lexmax()

        dimbox = dimmin.union(dimmax).bounded_simple_hull()
        fpdimmap = fpdimmap.intersect_range(dimbox)
        bbsym = bbsym.intersect(fpdimmap.domain())
    return bbsym
