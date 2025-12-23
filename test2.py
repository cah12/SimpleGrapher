import sympy as sp
from sympy import simplify, sqrt, sin, Piecewise, Interval, Expr
from sympy.solvers import solve


def find_discontinuities(expr, lower, upper):
    # Simplify the expression
    expr = simplify(expr)

    # Find the domain of the expression
    domain = solve(Expr(expr).as_real().eq(0), x)

    # Initialize the list of discontinuities
    discontinuities = []

    # Check each domain
    for d in domain:
        # Check if the domain is within the interval
        if Interval(lower, upper).contains(d):
            # Check the type of discontinuity
            if d in expr.free_symbols:
                # Essential discontinuity
                discontinuities.append({'position': d, 'type': 'essential'})
            else:
                # Check if the expression is defined at the domain
                if not expr.subs(x, d) in [0, 1]:
                    # Removable discontinuity
                    discontinuities.append(
                        {'position': d, 'type': 'removable', 'limit': expr.subs(x, d)})
                else:
                    # Check if the expression has a jump at the domain
                    if expr.limits(lower, upper, dir="+").evalf(subs={x: d}) != expr.limits(lower, upper, dir="-").evalf(subs={x: d}):
                        # Jump discontinuity
                        discontinuities.append({'position': d, 'type': 'jump', 'limit': expr.limits(
                            lower, upper, dir="+").evalf(subs={x: d}), 'limit2': expr.limits(lower, upper, dir="-").evalf(subs={x: d})})
                    else:
                        # Infinite discontinuity
                        discontinuities.append(
                            {'position': d, 'type': 'infinite'})

    return discontinuities


# quick example (uncomment to test)
if __name__ == "__main__":
    x = sp.symbols('x', real=True)
    expr = sp.sympify("sqrt(4+1/sin(x))")
    # expr = sp.sympify("1/sin(x)")
    print(find_discontinuities(expr,  -10, 10))
