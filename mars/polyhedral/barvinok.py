# This is dirty, but it's utterly complicated to get the Python interface to Barvinok
# (see https://repo.or.cz/barvinok.git). No easy way to get it as a PIP package seems to
# exist. Therefore we call Barvinok through a process (iscc) and grab the result into
# a piecewise quasi-polynomial (which is the most general kind of polynomial given by
# Barvinok)

import subprocess
from islpy import PwQPolynomial


def getBarvinokCardAsPWQPolynomial(s):
    if s.is_empty():
        return PwQPolynomial.read_from_str(s.get_ctx(), "{ 0 }")
    else:
        card_str = "card(" + s.to_str() + ");"
        result = subprocess.run(
            ["iscc"], input=card_str, capture_output=True, encoding="UTF-8"
        )
        assert result.returncode == 0, "Barvinok call failed"
        return PwQPolynomial.read_from_str(s.get_ctx(), result.stdout)


def getBarvinokCardAsInt(s):
    # Get a pwqpolynomial from barvinok...
    pwq = getBarvinokCardAsPWQPolynomial(s)
    pwq_pieces = pwq.get_pieces()
    assert (
        len(pwq_pieces) <= 1
    ), f"Cardinality of set is piecewise or not constant: {str(pwq)}"
    if len(pwq_pieces) == 0:
        return 0
    pwq_piece = pwq_pieces[0][1]
    return pwq_piece.get_constant_val().to_python()


def getBarvinokCardAsMinMaxInt(s):
    pwq = getBarvinokCardAsPWQPolynomial(s)
    pwq_pieces = pwq.get_pieces()
    M = max(
        [
            pwq_pieces[i][1].get_constant_val().to_python()
            for i in range(len(pwq_pieces))
        ]
    )
    m = min(
        [
            pwq_pieces[i][1].get_constant_val().to_python()
            for i in range(len(pwq_pieces))
        ]
    )
    assert m == M, f"Cardinality of set is not constant wrt parameters: {str(pwq)}"
    return m


def getBarvinokCardAsIntWithGist(s, param_rel):
    pwq = (
        getBarvinokCardAsPWQPolynomial(s.gist_params(param_rel.copy()))
        .gist_params(param_rel.copy())
        .coalesce()
    )
    pwq_pieces = pwq.get_pieces()
    M = max(
        [
            pwq_pieces[i][1].get_constant_val().to_python()
            for i in range(len(pwq_pieces))
        ]
    )
    m = min(
        [
            pwq_pieces[i][1].get_constant_val().to_python()
            for i in range(len(pwq_pieces))
        ]
    )
    assert m == M, f"Cardinality of set is not constant wrt parameters: {str(pwq)}"
    return m
