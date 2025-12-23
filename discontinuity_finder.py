import sympy as sp
from sympy import Symbol, Interval, Piecewise, EmptySet, limit, oo, nan, S, FiniteSet, Union, singularities


def _find_radical_discontinuities(
    expr: sp.Expr,
    var: sp.Symbol,
    lower: float,
    upper: float
) -> set:
    """Find discontinuities caused by even-order radicals."""
    radicals = set()

    # Find all radicals in the expression
    rad_terms = expr.find(sp.Pow)

    for term in rad_terms:
        # Check for even-order radicals or negative exponents
        base, exponent = term.as_base_exp()

        # Handle even roots like sqrt(f(x))
        if exponent.is_Rational:
            denom = exponent.q
            # Even root
            if denom % 2 == 0:
                # Find where the radicand is negative
                # radicand_zeros = solve(base, var)
                radicand_zeros = sp.solveset(base, var, sp.Interval(
                    lower, upper))
                for zero in radicand_zeros:
                    if zero.is_real:
                        point = float(zero)
                        if lower <= point <= upper:
                            radicals.add(zero)

    return radicals


def find_discontinuities(expr, x, lower, upper):
    """
    Finds all real discontinuities of expr in the interval [lower, upper] and their types.
    Returns a list of dictionaries: each with 'position' (float), 'type' (string), and optionally 'limit' (float) if exists.
    Types: 'removable', 'jump', 'infinite', 'essential'
    """
    search_interval = Interval(lower, upper)

    # 1. Find all singularities in the expression
    sing = singularities(expr, x)

    # 2. The discontinuities are the singularities within the search interval
    discontinuities = sing.intersection(search_interval)

    result = []

    def classify_point(c):
        try:
            left = limit(expr, x, c, '-')
            right = limit(expr, x, c, '+')
            f_c = expr.subs(x, c)
            if left == right and f_c == left:
                f_c1 = expr.subs(x, c+0.0001)
                f_c2 = expr.subs(x, c-0.0001)
                if f_c1.has(sp.I) and f_c2.has(sp.I):
                    return None, None
                return 'unknown2', f_c
            if left.has(sp.I) and right.has(sp.I):
                return 'unknown2', f_c
            if left == right:
                if left.is_finite:
                    if f_c == left:
                        return None, None  # continuous
                    else:
                        return 'removable', left
                else:
                    return 'infinite', None
            else:
                if left.is_finite and right.is_finite:
                    return 'jump', None
                else:
                    return 'essential', None
        except:
            return 'essential', None

    if discontinuities.is_FiniteSet:
        for c in discontinuities:
            typ, lim = classify_point(c)
            if typ:
                d = {'position': float(c), 'type': typ}
                if lim is not None:
                    d['limit'] = float(lim)
                result.append(d)
    elif discontinuities.is_Union:
        for part in discontinuities.args:
            if part.is_FiniteSet:
                for c in part:
                    typ, lim = classify_point(c)
                    if typ:
                        d = {'position': float(c), 'type': typ}
                        if lim is not None:
                            d['limit'] = float(lim)
                        result.append(d)

    # 3. Find discontinuities caused by radicals
    radicals = _find_radical_discontinuities(expr, x, lower, upper)
    for c in radicals:
        typ, lim = classify_point(c)
        if typ:
            d = {'position': float(c), 'type': typ}
            if lim is not None:
                d['limit'] = float(lim)
            result.append(d)

    # Skip intervals, as they represent undefined regions without specific points

    # Convert from a list of dictionaries to a list of lists
    result = [list(d.values()) for d in result]
    # Sort by position
    result.sort(key=lambda x: x[0])
    return result


# expr = sp.sympify("sqrt(sin(x))") #[]
# expr = sp.sympify("sqrt(1+sin(x))") #[]
# expr = sp.sympify("sqrt(4+sin(x))") #[]
# expr = sp.sympify("sqrt(4-sin(x))") #[]
# expr = sp.sympify("tan(x)") #[{'position': -7.853981633974483, 'type': 'essential'}, {'position': -4.71238898038469, 'type': 'essential'}, {'position': -1.5707963267948966, 'type': 'essential'}, {'position': 1.5707963267948966, 'type': 'essential'}, {'position': 4.71238898038469, 'type': 'essential'}, {'position': 7.853981633974483, 'type': 'essential'}]
# expr = sp.sympify("cot(x)") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.sympify("1/sin(x)") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.sympify("sqrt(1/sin(x))") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.sympify("sqrt(1+1/sin(x))") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
expr = sp.sympify("sqrt(4+1/sin(x))")  # [{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.sympify("abs(x)/x") #[{'position': 0.0, 'type': 'jump'}]
# expr = sp.sympify("sin(x)/x") #[{'position': 0.0, 'type': 'removable', 'limit': 1.0}]
# expr = sp.sympify("sqrt(1/x)") #[{'position': 0.0, 'type': 'removable', 'limit': 1.0}]
# expr = sp.sympify("(x+2)/(x^2-4)") #[[-2.0, 'removable', -0.25], [2.0, 'essential']]

# var = sp.sympify("x")

# discont = find_discontinuities(expr, var, -10, 10)

# print(discont)
# print(len(discont))
