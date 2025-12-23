import sympy as sp
from sympy import Interval, oo, S
from sympy.calculus.util import continuous_domain
from sympy import limit, singularities


def find_discontinuities(expr, x, lower, upper, tol=1e-8):
    """
    Finds real discontinuities of expr on [lower, upper].
    Returns a list of dicts: {'position': float, 'type': 'removable'|'jump'|'infinite'|'essential', 'limit': float (optional)}
    Designed to handle nested domain issues (e.g. sqrt(4+1/sin(x))).
    """
    search = Interval(lower, upper)

    # domain where expression is real and continuous on the search interval
    try:
        domain = continuous_domain(expr, x, search)
    except Exception:
        domain = search  # conservative fallback

    # the excluded set inside the search interval
    excluded = search - domain

    # include symbolic singularities (e.g. poles) as candidates
    try:
        sing = singularities(expr, x)
    except Exception:
        sing = S.EmptySet

    candidates = set()

    def add_sym_point(p):
        if p is None or p is S.EmptySet:
            return
        # only real, finite candidates inside the search interval
        try:
            if p.is_real and p.is_finite and p in search:
                candidates.add(sp.nsimplify(p))
        except Exception:
            # nsimplify may fail for transcendental points; try numeric check
            try:
                pv = float(sp.N(p))
                if lower - tol <= pv <= upper + tol:
                    candidates.add(sp.N(p))
            except Exception:
                pass

    # gather from excluded (FiniteSet or Interval endpoints)
    if excluded is not S.EmptySet:
        if excluded.is_FiniteSet:
            for p in excluded:
                add_sym_point(p)
        elif excluded.is_Union:
            for part in excluded.args:
                if part.is_FiniteSet:
                    for p in part:
                        add_sym_point(p)
                elif part.is_Interval:
                    for endpoint in (part.start, part.end):
                        add_sym_point(endpoint)
        elif excluded.is_Interval:
            for endpoint in (excluded.start, excluded.end):
                add_sym_point(endpoint)

    # gather from singularities
    if sing is not S.EmptySet:
        if getattr(sing, "is_FiniteSet", False):
            for p in sing:
                if p in search:
                    add_sym_point(p)
        elif getattr(sing, "is_Union", False):
            for part in sing.args:
                if getattr(part, "is_FiniteSet", False):
                    for p in part:
                        if p in search:
                            add_sym_point(p)

    # helper to test symbolic infinities
    def is_symbolic_infinite(v):
        return v in (sp.oo, -sp.oo, sp.zoo, S.ComplexInfinity)

    # numeric extractor: returns (float_value, bad_flag, is_infinite_flag, is_complex_flag)
    def numeric_value(v):
        try:
            if is_symbolic_infinite(v):
                return None, True, True, False
            vn = sp.N(v, 20)
            # complex results
            if getattr(vn, "as_real_imag", None):
                re, im = sp.re(vn), sp.im(vn)
                try:
                    imf = float(im.evalf())
                except Exception:
                    return None, True, False, True
                if abs(imf) > tol:
                    return None, True, False, True
                rf = float(re.evalf())
                if abs(rf) > 1e300:
                    return None, True, True, False
                return rf, False, False, False
            # fallback
            rf = float(vn)
            if abs(rf) > 1e300:
                return None, True, True, False
            return rf, False, False, False
        except Exception:
            return None, True, False, True

    result = []

    for p in sorted(candidates, key=lambda z: float(sp.N(z))):
        # use exact sympy point for limits where possible
        c = p
        try:
            L = limit(expr, x, c, dir='-')
            R = limit(expr, x, c, dir='+')
        except Exception:
            # if symbolic limit fails, mark essential
            result.append({'position': float(sp.N(c)), 'type': 'essential'})
            continue

        # quick check for symbolic infinities
        if is_symbolic_infinite(L) or is_symbolic_infinite(R):
            result.append({'position': float(sp.N(c)), 'type': 'infinite'})
            continue

        # numeric analysis
        lval, lbad, linf, lcomp = numeric_value(L)
        rval, rbad, rinf, rcomp = numeric_value(R)

        # if either side is infinite symbolically or numerically
        if linf or rinf:
            result.append({'position': float(sp.N(c)), 'type': 'infinite'})
            continue

        # if either limit is complex or cannot be numeric -> essential
        if lcomp or rcomp or lbad or rbad:
            result.append({'position': float(sp.N(c)), 'type': 'essential'})
            continue

        # both sides numeric finite
        if abs(lval - rval) <= max(tol, tol * max(abs(lval), abs(rval), 1.0)):
            # check function value at point
            try:
                fc = expr.subs(x, c)
                fnum, fbad, finf, fcomp = numeric_value(fc)
                if not fbad and not fcomp and not finf and abs(fnum - lval) <= tol:
                    # continuous
                    continue
                else:
                    # removable (limits equal but function undefined or different)
                    entry = {'position': float(
                        sp.N(c)), 'type': 'removable', 'limit': float(lval)}
                    result.append(entry)
                    continue
            except Exception:
                entry = {'position': float(
                    sp.N(c)), 'type': 'removable', 'limit': float(lval)}
                result.append(entry)
                continue
        else:
            # finite but different one-sided limits -> jump
            result.append({'position': float(sp.N(c)), 'type': 'jump'})
            continue

    # sort by position
    result.sort(key=lambda d: d['position'])
    return result


# quick example (uncomment to test)
if __name__ == "__main__":
    x = sp.symbols('x', real=True)
    expr = sp.sympify("sqrt(4+1/sin(x))")
    # expr = sp.sympify("1/sin(x)")
    print(find_discontinuities(expr, x, -10, 10))
