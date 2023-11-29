# Warning: the lists that represent sets have to be sorted the same way
# (list comparison matches the order while set doesn't)


def partsOfSet(set_lst):
    def diff(a, b):
        return [x for x in a if not x in b]

    def add(a, b):
        return a + b

    def parts_rec(set_arg):
        if len(set_arg) == 0:
            return [[]]
        for a in set_arg:
            subsets = parts_rec(diff(set_arg, [a]))
            return list(map(lambda s: s + [a], subsets)) + subsets

    return parts_rec(set_lst)


def nonTrivialPartsOfSet(set_lst):
    return list(filter(lambda l: l != [], partsOfSet(set_lst)))


def scalarProduct(a, b):
    assert len(a) == len(b)
    s = 0
    for i in range(len(a)):
        s += a[i] * b[i]
    return s


def complement(A, B):
    return [x for x in B if not x in A]


def listEqualsUnordered(A, B):
    if len(A) != len(B):
        return False
    for a in A:
        if a not in B:
            return False
    return True
