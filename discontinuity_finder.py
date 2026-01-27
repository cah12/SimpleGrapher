import sympy as sp
from sympy import lambdify, Symbol, Interval, Piecewise, EmptySet, limit, oo, nan, S, FiniteSet, Union, singularities
from degree_radian import get_mode
from solveset_thread import limit_with_timeout
import numpy as np
from sympy.functions.elementary.trigonometric import TrigonometricFunction


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
        if term.is_real:
            continue
        # Check for even-order radicals or negative exponents
        base, exponent = term.as_base_exp()

        # Handle even roots like sqrt(f(x))
        if exponent.is_Rational or exponent.is_rational:
            if exponent.is_Rational:
                denom = exponent.q
            else:
                num, denom = exponent.as_numer_denom()
            # Even root
            if denom % 2 == 0:
                # Find where the radicand is negative
                # radicand_zeros = solve(base, var)
                radicand_zeros = sp.solveset(base, var, sp.Interval(
                    lower, upper))
                try:
                    for zero in radicand_zeros:
                        if zero.is_real:
                            point = float(zero)
                            if lower <= point <= upper:
                                radicals.add(zero)
                except:
                    pass

    return radicals


def nsolveEquation(equation, start, stop, var, numOfGuess=300, decimalPlaces=4):
    # 2. Define the range for initial guesses using linspace
    # We need to ensure the guesses are within the domain where roots exist
    start, stop = start, stop
    guesses = np.linspace(start, stop, numOfGuess)

    # 3. Use nsolve in a loop to find roots for each initial guess
    roots = set()

    for guess in guesses:
        try:
            root = sp.nsolve(equation, var, guess,
                             verify=True, prec=decimalPlaces)
            if sp.im(root) == 0:  # Only consider real roots
                root = sp.re(root)
                if float(root) >= start and float(root) <= stop:
                    roots.add(root)
                    continue
        except (ValueError, RuntimeError):
            # nsolve may fail if the initial guess is poor or a singularity is hit
            continue

    # Filter unique roots if necessary
    return roots


def find_discontinuities(expr, x, lower, upper, period):
    """
    Finds all real discontinuities of expr in the interval [lower, upper] and their types.
    Returns a list of dictionaries: each with 'position' (float), 'type' (string), and optionally 'limit' (float) if exists.
    Types: 'removable', 'jump', 'infinite', 'essential'
    """

    def classify_point(c):
        try:
            # if not period and c.is_real:
            if not expr.has(TrigonometricFunction) and c.is_real:
                c = float(c)
            left = limit(expr, x, c, '-')
            right = limit(expr, x, c, '+')

            f_c = expr.subs(x, c)
            f_c = sp.simplify(f_c)

            if not expr.has(TrigonometricFunction):
                left = left.evalf()
                right = right.evalf()
                f_c = f_c.evalf()

            if f_c.is_real:
                f_c = f_c

            if left == right and f_c == left:
                f_c1 = expr.subs(x, c+0.0001)
                f_c2 = expr.subs(x, c-0.0001)
                if f_c1.has(sp.I) and f_c2.has(sp.I):
                    return None, None
                return 'unknown2', f_c
            if left.has(sp.I) and right.has(sp.I):
                if not f_c.is_real:
                    f_c = None
                return 'unknown2', f_c
            if left == right:
                if left.is_finite:
                    if not left.args and f_c == float(left):
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

    result = []

    search_interval = Interval(lower, upper)

    # 1. Find all singularities in the expression
    sing = singularities(expr, x)

    # 2. The discontinuities are the singularities within the search interval
    discontinuities = sing.intersection(search_interval)

    if isinstance(discontinuities, sp.ConditionSet):
        # If singularities fails, try to find them manually
        try:
            # For rational functions, find zeros of denominator
            numer, denom = sp.fraction(expr)
            # denom_zeros = sp.solve(denom, x)
            discontinuities = sp.solveset(sp.simplify(denom), x, sp.Interval(
                lower, upper))  # solve(denom, var)
            for part in discontinuities.args:
                if isinstance(part, sp.Intersection):
                    list_of_lambda = [
                        r for r in part.args if isinstance(r, sp.ImageSet)]
                    _n = sp.symbols("_n")
                    for c in list_of_lambda:
                        ld = c.args[0]
                        ld_expr = str(ld.expr)
                        ld_expr = sp.parse_expr(ld_expr)

                        for n in range(0, 50):
                            ld_expr_subs = ld_expr.subs(_n, n)
                            if not ld_expr_subs.is_real:
                                break
                            if float(ld_expr_subs) <= lower or float(ld_expr_subs) >= upper:
                                break

                            d = {'position': float(
                                ld_expr_subs), 'type': "essential"}
                            result.append(d)

                if isinstance(part, sp.Eq):
                    e = part.args[0]
                    dis = nsolveEquation(e, -10, 10, x)

                    for c in dis:
                        if not c.is_real:
                            continue
                        if float(c) <= lower or float(c) >= upper:
                            break
                        d = {'position': float(c), 'type': "essential"}
                        result.append(d)
                    break
        except Exception:
            pass

    # 3. Classify the discontinuities

    # result = []

    if discontinuities.is_FiniteSet:
        for c in discontinuities:
            typ, lim = classify_point(c)
            if typ:
                d = {'position': float(c), 'type': typ}
                if lim is not None:
                    try:
                        lim = float(lim)
                        d['limit'] = float(lim)
                    except:
                        pass
                result.append(d)
    elif discontinuities.is_Union:
        # 1.2533141,2.8024956
        for part in discontinuities.args:
            if part.is_FiniteSet:
                for c in part:
                    if not c.is_real:
                        continue
                    if float(c) <= lower or float(c) >= upper:
                        break
                    typ, lim = classify_point(c)
                    if typ:
                        d = {'position': float(c), 'type': typ}
                        if lim is not None:
                            d['limit'] = float(lim)
                        result.append(d)
            elif isinstance(part, sp.ConditionSet):
                eq = part.args[1]
                dis = nsolveEquation(eq, lower, upper, x)
                for c in dis:
                    if not c.is_real:
                        continue
                    # typ, lim = classify_point(c)
                    if float(c) <= lower or float(c) >= upper:
                        break
                    d = {'position': float(c), 'type': "essential"}
                    result.append(d)

    # 3. Find discontinuities caused by radicals
    radicals = _find_radical_discontinuities(expr, x, lower, upper)
    for c in radicals:
        typ, lim = classify_point(c)
        if typ:
            d = {'position': float(c), 'type': typ}
            if lim is not None:
                try:
                    d['limit'] = float(lim)
                except:
                    if not lim.is_real and lim.args[0].is_real:
                        d['limit'] = float(lim.args[0])

            result.append(d)

    # Skip intervals, as they represent undefined regions without specific points

    # Convert from a list of dictionaries to a list of lists
    result = [list(d.values()) for d in result]

    # Sort by position
    result.sort(key=lambda x: x[0])
    return result


# expr = sp.parse_expr("sqrt(sin(x))") #[]
# expr = sp.parse_expr("sqrt(1+sin(x))") #[]
# expr = sp.parse_expr("sqrt(4+sin(x))") #[]
# expr = sp.parse_expr("sqrt(4-sin(x))") #[]
# expr = sp.parse_expr("tan(x)") #[{'position': -7.853981633974483, 'type': 'essential'}, {'position': -4.71238898038469, 'type': 'essential'}, {'position': -1.5707963267948966, 'type': 'essential'}, {'position': 1.5707963267948966, 'type': 'essential'}, {'position': 4.71238898038469, 'type': 'essential'}, {'position': 7.853981633974483, 'type': 'essential'}]
# expr = sp.parse_expr("cot(x)") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.parse_expr("1/sin(x)") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.parse_expr("sqrt(1/sin(x))") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.parse_expr("sqrt(1+1/sin(x))") #[{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# [{'position': 0.0, 'type': 'essential'}, {'position': 3.141592653589793, 'type': 'essential'}, {'position': -9.42477796076938, 'type': 'essential'}, {'position': -6.283185307179586, 'type': 'essential'}, {'position': -3.141592653589793, 'type': 'essential'}, {'position': 6.283185307179586, 'type': 'essential'}, {'position': 9.42477796076938, 'type': 'essential'}]
# expr = sp.parse_expr("sqrt(4+1/sin(x))")
# expr = sp.parse_expr("abs(x)/x") #[{'position': 0.0, 'type': 'jump'}]
# expr = sp.parse_expr("sin(x)/x") #[{'position': 0.0, 'type': 'removable', 'limit': 1.0}]
# expr = sp.parse_expr("sqrt(1/x)") #[{'position': 0.0, 'type': 'removable', 'limit': 1.0}]
# expr = sp.parse_expr("(x+2)/(x^2-4)") #[[-2.0, 'removable', -0.25], [2.0, 'essential']]

# var = sp.parse_expr("x")

# discont = find_discontinuities(expr, var, -10, 10)

# print(discont)
# print(len(discont))
