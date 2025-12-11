"""
Module for finding discontinuities in SymPy expressions.
"""

import sympy as sp
from sympy import symbols, limit, oo, re, im, diff, solve, simplify
# from sympy.analysis.singularities import singularities
from sympy.calculus.singularities import singularities
from typing import List, Tuple, Dict, Optional

#  This code defines a function
# find_discontinuities
#  that takes in an expression (expr), a variable (var), and a range of values (lower_limit and upper_limit). The function identifies and classifies discontinuities in the expression within the specified range. The discontinuities are classified as "infinite", "removable", "jump", or "unknown". The function handles various expression types including rational functions, trigonometric functions, radicals, and logarithms.

# The function first checks if the lower limit is less than the upper limit. If not, it raises a ValueError. It then finds candidate points where discontinuities may occur using helper functions _find_pole_points,
# _find_radical_discontinuities
# ,
# _find_logarithm_discontinuities
# , and
# _find_trig_discontinuities
# c:\Users\chope\Downloads\grapher-master\grapher-master\discontinuity_finder.py
# . It then classifies each discontinuity using the
# _classify_discontinuity
#  function and appends the classification to a list. Finally, it returns the list of discontinuity classifications.

# The function uses the SymPy library to perform the calculations and analysis. The return value is a list of dictionaries, where each dictionary contains information about a discontinuity, such as the point where it occurs, the type of discontinuity, and the left and right limits (if applicable).


def find_discontinuities(
    expr: sp.Expr,
    var: sp.Symbol,
    lower_limit: float,
    upper_limit: float,
) -> List[Dict[str, any]]:
    """
    Find all discontinuities of a SymPy expression within a specified range.

    This function identifies and classifies discontinuities as:
    - "infinite": Vertical asymptotes where the function diverges
    - "removable": Points where the limit exists but the function is undefined
    - "jump": Points where left and right limits exist but differ
    - "unknown": Points where discontinuity type cannot be determined

    The function handles various expression types including rational functions,
    trigonometric functions, radicals, and logarithms.

    Parameters
    ----------
    expr : sympy.Expr
        The expression to analyze for discontinuities
    var : sympy.Symbol
        The variable with respect to which to find discontinuities
    lower_limit : float
        Lower bound of the range to search (inclusive)
    upper_limit : float
        Upper bound of the range to search (inclusive)

    Returns
    -------
    List[Dict[str, any]]
        A list of dictionaries, each containing:
        - "point": The x-coordinate where discontinuity occurs (float)
        - "type": Classification as "infinite", "removable", "jump", or "unknown"
        - "left_limit": Left-hand limit (or None if infinite)
        - "right_limit": Right-hand limit (or None if infinite)
        - "limit": The limit value for removable discontinuities (or None otherwise)

    Examples
    --------
    >>> from sympy import symbols, sin, pi
    >>> x = symbols('x')
    >>> 
    >>> # Find discontinuities of 1/x between -5 and 5
    >>> discontinuities = find_discontinuities(1/x, x, -5, 5)
    >>> # Returns: [{'point': 0.0, 'type': 'infinite', ...}]
    >>> 
    >>> # Find discontinuities of tan(x) between -π and π
    >>> discontinuities = find_discontinuities(sp.tan(x), x, -sp.pi, sp.pi)
    >>> # Returns discontinuities at x = ±π/2
    """

    if isinstance(expr, str):
        expr = sp.sympify(expr)

    # Convert limits to floats for comparison
    lower = float(lower_limit)
    upper = float(upper_limit)

    if lower >= upper:
        raise ValueError("lower_limit must be less than upper_limit")

    discontinuities = []
    candidate_points = set()

    # Find candidate discontinuity points
    candidate_points.update(_find_pole_points(expr, var, lower, upper))
    candidate_points.update(
        _find_radical_discontinuities(expr, var, lower, upper))
    candidate_points.update(
        _find_logarithm_discontinuities(expr, var, lower, upper))
    candidate_points.update(
        _find_trig_discontinuities(expr, var, lower, upper))

    # Classify each discontinuity
    for point in sorted(candidate_points):
        if lower <= point <= upper:
            classification = _classify_discontinuity(expr, var, point)
            if classification:
                discontinuities.append(classification)

    return discontinuities


""" This code defines a private function 
_find_pole_points
 that takes in an expression, a symbol, and a range of values. It returns a set of points where the expression has infinite discontinuities, also known as poles.

The function first initializes an empty set called poles to store the points of discontinuity. It then tries to find the singularities of the expression using the singularities method provided by the SymPy library. It iterates over each singularity and checks if it is a real number within the given range. If it is, the point is added to the poles set.

If the singularities method fails, the function tries to find the zeros of the denominator of the expression using the fraction method provided by SymPy. It then iterates over each zero and checks if it is a real number within the range. If it is, the point is added to the poles set.

Finally, the function retur """


def _find_pole_points(
    expr: sp.Expr,
    var: sp.Symbol,
    lower: float,
    upper: float
) -> set:
    """Find poles (infinite discontinuities) from singularities."""
    poles = set()

    try:
        # Get singularities of the expression
        sings = singularities(expr, var)

        for sing in sings:
            # Check if it's a real number within range
            if sing.is_real:
                point = float(sing)
                if lower <= point <= upper:
                    poles.add(point)
    except Exception:
        # If singularities fails, try to find them manually
        try:
            # For rational functions, find zeros of denominator
            numer, denom = sp.fraction(expr)
            denom_zeros = solve(denom, var)

            for zero in denom_zeros:
                if zero.is_real:
                    point = float(zero)
                    if lower <= point <= upper:
                        poles.add(point)
        except Exception:
            pass

    return poles


""" This code snippet defines a private function 
_find_radical_discontinuities
 that takes in an expression, a symbol, and a range of values. It returns a set of points where the logarithmic terms in the expression are undefined within the given range.

The function first initializes an empty set called 
radical_discontinuities
 to store the points of discontinuity. It then uses the 
find
 method provided by the SymPy library to find all occurrences of the Pow function in the expression.

Next, it iterates over each term found and retrieves the base and exponent of the term using the as_base_exp method. It checks if the exponent is a rational number.

If the exponent is a rational number, it checks if the denominator of the rational number is even. If so, it uses the solve method to find the values of the variable where the base equals zero.

For each zero found, it checks if it is a real number and if it falls within the specified range. If so, it adds the point to the 
radical_discontinuities
 set.

Finally, it returns the set of points where the logarithmic terms are undefined within the given range. """


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
                radicand_zeros = solve(base, var)
                for zero in radicand_zeros:
                    if zero.is_real:
                        point = float(zero)
                        if lower <= point <= upper:
                            radicals.add(point)

    return radicals


""" This code defines a private function 
_find_logarithm_discontinuities
 that takes in an expression, a symbol, and a range of values. It returns a set of points where the logarithmic terms in the expression are undefined within the given range.

The function initializes an empty set called logs to store the points of discontinuity. It then uses the 
find
 method provided by the SymPy library to find all occurrences of the logarithm function (sp.log) in the expression.

Next, it iterates over each term found and retrieves the argument of the logarithm function. It then uses the solve method to find the values of the variable where the argument equals zero.

For each zero found, it checks if it is a real number and if it falls within the specified range. If so, it adds the point to the logs set.

Finally, it returns the set of points of discontinuity. """


def _find_logarithm_discontinuities(
    expr: sp.Expr,
    var: sp.Symbol,
    lower: float,
    upper: float
) -> set:
    """Find discontinuities caused by logarithmic terms."""
    logs = set()

    # Find all logarithms in the expression
    log_terms = expr.find(sp.log)

    for term in log_terms:
        arg = term.args[0]
        # Find where the argument equals zero
        zeros = solve(arg, var)

        for zero in zeros:
            if zero.is_real:
                point = float(zero)
                if lower <= point <= upper:
                    logs.add(point)

    return logs


""" This code snippet defines a private function 
_find_trig_discontinuities
 that finds discontinuities caused by trigonometric functions (tan, cot, sec, csc).

The function takes in an expression (expr), a symbol (var), and a range of values (lower and upper). It returns a set of points where the trigonometric functions are undefined within the given range.

The function first initializes an empty set trig_discs to store the points of discontinuity. It then searches for occurrences of tan, cot, sec, and csc functions in the expression using the 
find
 method provided by SymPy library.

Next, it iterates over each term found and checks if it is either tan or sec. If so, it solves for the points where the argument (arg) of the term is equal to π/2 + nπ for n ranging from -10 to 10. It then checks if the solution is real and within the given range. If so, it adds the point to the trig_discs set.

If the term is cot or csc, it solves for the points where the argument (arg) of the term is equal to nπ for n ranging from -10 to 10. It then checks if the solution is real and within the given range. If so, it adds the point to the trig_discs set.

Finally, it returns the set of points of discontinuity. """


def _find_trig_discontinuities(
    expr: sp.Expr,
    var: sp.Symbol,
    lower: float,
    upper: float
) -> set:
    """Find discontinuities caused by trigonometric functions (tan, cot, sec, csc)."""
    trig_discs = set()

    # Check for tan, cot, sec, csc
    tan_terms = expr.find(sp.tan)
    cot_terms = expr.find(sp.cot)
    sec_terms = expr.find(sp.sec)
    csc_terms = expr.find(sp.csc)

    all_trig = tan_terms | cot_terms | sec_terms | csc_terms

    for term in all_trig:
        arg = term.args[0]

        # Find discontinuities based on function type
        if isinstance(term, sp.tan) or isinstance(term, sp.sec):
            # tan(x) and sec(x) undefined at π/2 + nπ
            # Solve: arg = π/2 + nπ
            for n in range(-10, 11):
                discontinuity_point = sp.pi / 2 + n * sp.pi
                # Solve arg = discontinuity_point
                solutions = solve(arg - discontinuity_point, var)
                for sol in solutions:
                    if sol.is_real:
                        point = float(sol)
                        if lower <= point <= upper:
                            trig_discs.add(point)

        elif isinstance(term, sp.cot) or isinstance(term, sp.csc):
            # cot(x) and csc(x) undefined at nπ
            for n in range(-10, 11):
                discontinuity_point = n * sp.pi
                solutions = solve(arg - discontinuity_point, var)
                for sol in solutions:
                    if sol.is_real:
                        point = float(sol)
                        if lower <= point <= upper:
                            trig_discs.add(point)

    return trig_discs


""" This code snippet defines a function 
_classify_discontinuity
 that takes an expression (expr), a symbol (var), and a point (point) as input. It classifies a discontinuity at the given point based on the behavior of the expression.

The function first evaluates the left and right limits of the expression at the given point using the limit function from the SymPy library. It then tries to evaluate the expression at the point and checks if it is defined and finite.

Based on the values of the left and right limits, the function classifies the discontinuity into different types:

If either the left or right limit is infinite, it is an infinite discontinuity.
If the left and right limits are equal and the expression is not defined at the point, it is a removable discontinuity.
If the left and right limits are different and neither is infinite, it is a jump discontinuity.
If the left and right limits are different or the expression is not equal to the point value, it is an unknown or other type of discontinuity.
The function returns a dictionary containing the classification details, including the point, type of discontinuity, left limit, right limit, and the limit value (if applicable). If no discontinuity is found, it returns None. """


def _classify_discontinuity(
    expr: sp.Expr,
    var: sp.Symbol,
    point: float
) -> Optional[Dict[str, any]]:
    """
    Classify a discontinuity at a given point.

    Returns a dictionary with classification details or None if no discontinuity.
    """

    # Evaluate left and right limits
    left_lim = limit(expr, var, point, '-')
    right_lim = limit(expr, var, point, '+')

    # Try to evaluate the function at the point
    try:
        point_value = expr.subs(var, point)
        # Check if it's defined and finite
        is_defined = not (point_value.has(sp.zoo) or point_value.has(sp.oo) or
                          point_value.has(-sp.oo) or point_value is sp.nan)
    except Exception:
        is_defined = False

    # Classify the discontinuity
    if left_lim == oo or left_lim == -oo or right_lim == oo or right_lim == -oo:
        # Infinite discontinuity
        return {
            "point": point,
            "type": "infinite",
            "left_limit": None if left_lim in (oo, -oo) else float(left_lim) if left_lim.is_number else left_lim,
            "right_limit": None if right_lim in (oo, -oo) else float(right_lim) if right_lim.is_number else right_lim,
            "limit": None
        }

    elif left_lim == right_lim and not is_defined:
        # Removable discontinuity
        if left_lim.is_number:
            limit_value = float(left_lim)
        else:
            limit_value = None

        return {
            "point": point,
            "type": "removable",
            "left_limit": limit_value,
            "right_limit": limit_value,
            "limit": limit_value
        }

    elif left_lim != right_lim and left_lim not in (oo, -oo) and right_lim not in (oo, -oo):
        # Jump discontinuity
        left_val = float(left_lim) if left_lim.is_number else left_lim
        right_val = float(right_lim) if right_lim.is_number else right_lim

        return {
            "point": point,
            "type": "jump",
            "left_limit": left_val,
            "right_limit": right_val,
            "limit": None
        }

    elif left_lim != right_lim or (is_defined and left_lim != point_value):
        # Unknown or other type of discontinuity
        return {
            "point": point,
            "type": "unknown",
            "left_limit": float(left_lim) if left_lim.is_number else None,
            "right_limit": float(right_lim) if right_lim.is_number else None,
            "limit": None
        }

    return None


# Example usage and testing
if __name__ == "__main__":
    x = sp.symbols('x')

    print("=" * 70)
    print("Example 1: f(x) = 1/x")
    print("=" * 70)
    expr1 = 1 / x
    discs1 = find_discontinuities(expr1, x, -5, 5)
    for disc in discs1:
        print(f"Point: {disc['point']}, Type: {disc['type']}")
        if disc['type'] == 'removable':
            print(f"  Limit: {disc['limit']}")

    print("\n" + "=" * 70)
    print("Example 2: f(x) = sin(x)/x")
    print("=" * 70)
    expr2 = sp.sin(x) / x
    discs2 = find_discontinuities(expr2, x, -5, 5)
    for disc in discs2:
        print(f"Point: {disc['point']}, Type: {disc['type']}")
        if disc['type'] == 'removable':
            print(f"  Limit: {disc['limit']}")

    print("\n" + "=" * 70)
    print("Example 3: f(x) = tan(x)")
    print("=" * 70)
    expr3 = sp.tan(x)
    discs3 = find_discontinuities(expr3, x, -sp.pi, sp.pi)
    for disc in discs3:
        print(f"Point: {disc['point']:.4f}, Type: {disc['type']}")

    print("\n" + "=" * 70)
    print("Example 4: f(x) = sqrt(x-1)")
    print("=" * 70)
    expr4 = sp.sqrt(x - 1)
    discs4 = find_discontinuities(expr4, x, -2, 5)
    for disc in discs4:
        print(f"Point: {disc['point']}, Type: {disc['type']}")

    print("\n" + "=" * 70)
    print("Example 5: f(x) = (x^2 - 1)/(x - 1)")
    print("=" * 70)
    expr5 = (x**2 - 1) / (x - 1)
    discs5 = find_discontinuities(expr5, x, -3, 3)
    for disc in discs5:
        print(f"Point: {disc['point']}, Type: {disc['type']}")
        if disc['type'] == 'removable':
            print(f"  Limit: {disc['limit']}")
