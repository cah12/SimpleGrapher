from ast import expr
from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

import sympy as sp
from sympy import symbols, solve, fraction, oo, S, preorder_traversal, sin, cos, tan, cot, sec, csc
from sympy.functions.elementary.trigonometric import TrigonometricFunction
from sympy.core.function import _mexpand as flat
from typing import List, Union, Tuple
from sympy.calculus.util import periodicity

from degree_radian import trig_substitutions, set_mode, get_mode
from discontinuity_finder import find_discontinuities

from solveset_thread import solve_with_timeout


# The Degree To Radian Code
# --------------------------#
####################################################
###########################################################

sympified = False
full_expr = ""

""" expr = expr.subs(sp.sin, sin_mode)
    expr = expr.subs(sp.cos, cos_mode) """


def pyExpToJsExp(s):
    while s.find("**") != -1:
        s = s.replace("**", "^")
    return s


def jsExpToPyExp(s):
    while s.find("^") != -1:
        s = s.replace("^", "**")
    return s


def solve_for(exp, c):
    # while exp.find("^") != -1:
    #     exp = exp.replace("^", "**")

    exp = jsExpToPyExp(exp)

    v = sp.Symbol(c)
    arr = exp.split('=')

    if len(arr) == 2:
        if arr[0] == c and c not in arr[1]:
            return [pyExpToJsExp(arr[1])]
        if arr[1] == c and c not in arr[0]:
            return [pyExpToJsExp(arr[0])]

    if len(arr) == 1:
        arr.append("0")

    result = []
    try:
        solutions = None
        arr0 = None
        arr1 = None
        # with sp.evaluate(False):
        arr0 = sp.sympify(arr[0])
        if arr[1] == "0":
            solutions = sp.solveset(arr0, v)
        else:
            # with sp.evaluate(False):
            arr1 = sp.sympify(arr[1])
            s = sp.solveset(
                sp.Eq(arr0, arr1), v)
            if s.is_FiniteSet:
                solutions = list(s)
            elif s.args:
                for s in s.args:
                    if s.is_FiniteSet:
                        solutions = list(s)
                        break

        l = len(solutions)

        if l > 2:
            if l % 2 == 0:
                solutions = [solutions[0], solutions[1]]
            else:
                solutions = [solutions[0]]

        for solution in solutions:
            num, denom = solution.as_numer_denom()
            if denom != 1:
                if num.is_polynomial():
                    num = sp.factor(num)

                if denom.is_polynomial():
                    denom = sp.factor(denom)

                solution = (num)/(denom)
                solution = solution.simplify()

            _str = str(solution)
            # while _str.find("**") != -1:
            #     _str = _str.replace("**", "^")
            _str = pyExpToJsExp(_str)
            if _str.find("Piecewise") == -1 and _str.find("I") == -1:
                result.append(_str)

    except BaseException as error:
        print(error)
        return []

    return result


def inflection_points(expr, lower, upper, var):
    return []


def turning_points2(expr, lower, upper, var):
    try:
        d = sp.diff(sp.parse_expr(expr), sp.Symbol(var))
        dis = solve_with_timeout(d, var, sp.Interval(
            lower, upper, left_open=True, right_open=True))
        return list(map(float, dis))
    except:
        return []


def turning_points(expression, dependent_variable, lower_limit, upper_limit):
    var = sp.Symbol(dependent_variable)

    parsed_expression = sp.sympify(expression)

    # parsed_expression = trig_substitutions(parsed_expression)
    # if get_mode() == "deg":
    #     lower_limit = lower_limit * sp.pi / 180
    #     upper_limit = upper_limit * sp.pi / 180

    f_prime = sp.diff(parsed_expression, var)
    f_prime = trig_substitutions(f_prime)
    parsed_expression = trig_substitutions(parsed_expression)

    if isinstance(f_prime, sp.Number):
        return []

    points = []

    try:
        critical_points = solve_with_timeout(f_prime, var, sp.Interval(
            lower_limit, upper_limit, left_open=True, right_open=True))

        # critical_points = sp.solveset(f_prime, var, sp.Interval(
        #     lower_limit, upper_limit, left_open=True, right_open=True))

        num = 0
        if not isinstance(critical_points, sp.Complement):
            for point in critical_points:
                if num > 200:
                    break
                if point < lower_limit or point > upper_limit:
                    num += 1
                    continue
                max_x = point
                if max_x.is_real:
                    max_y = parsed_expression.subs(var, max_x)
                if max_y.is_real:
                    points.append([float(max_x), float(max_y)])

                num += 1
    except:
        pass

    return points


def inflection_points(expr, lower, upper, var):
    try:
        x = sp.Symbol(var)
        d = sp.diff(sp.parse_expr(expr), x, x)
        dis = solve_with_timeout(d, var, sp.Interval(lower, upper))
        return list(map(float, dis))
    except:
        return []


def handlePeriodic(discontinuitiesArr, lower, upper):
    if not isinstance(discontinuitiesArr, list) or len(discontinuitiesArr) < 2:
        return discontinuitiesArr
    result = []

    d = 2*sp.pi

    for discont in discontinuitiesArr:
        if analyze_discontinuity_type(full_expr, discont) == "unknown":
            continue
        if not discont.is_real:
            continue
        a1 = discont
        result.append(a1)
        a1 = a1 - d
        while a1 > lower:
            result.append(a1)
            a1 = a1 - d

        a1 = discont
        a1 = a1 + d
        while a1 < upper:
            result.append(a1)
            a1 = a1 + d

    result.sort()
    result = unique_elements(result)
    return result


""" static handlePeriodic(discontinuitiesArr, lower, upper) {
    if (!Array.isArray(discontinuitiesArr) | | discontinuitiesArr.length < 2) {
        return discontinuitiesArr
    }
    let d = discontinuitiesArr[1][0] - discontinuitiesArr[0][0]
    let a1 = discontinuitiesArr[0][0]

    if (d != 0) {

        let n = 0
        while (a1 > lower & & n < 5000) {
            a1 = a1 - n * d
            n++
        }
        n = 0
        while (a1 < lower & & n < 5000) {
            a1 = a1 + n * d
            n++
        }

        a1 = a1 - 3 * d
        n = 0
        discontinuitiesArr.length = 0

        while (n < 5000) {
            discontinuitiesArr.push([a1 + n * d, "infinite"])
            if (discontinuitiesArr[n][0] > upper) {
                break
            }
            n++
        }
        discontinuitiesArr.push([a1 + (n + 1) * d, "infinite"])
        discontinuitiesArr.push([a1 + (n + 2) * d, "infinite"])
        discontinuitiesArr.push([a1 + (n + 3) * d, "infinite"])
    }
    return discontinuitiesArr
} """


def find_discontinuities_in_range(
    expr: Union[str, sp.Expr],
    x_min: float,
    x_max: float,
    var: Union[str, sp.Symbol] = None
) -> List[float]:
    """
    Find all discontinuities in an algebraic expression within a given range [x_min, x_max].

    Discontinuities occur when:
    1. The denominator equals zero (division by zero)
    2. Logarithms have non-positive arguments
    3. Even roots have negative arguments
    4. Other undefined points

    Parameters:
    -----------
    expr : str or sympy expression
        The algebraic expression to analyze
    x_min : float
        Lower bound of the range (inclusive)
    x_max : float
        Upper bound of the range (inclusive)
    var : str or sympy Symbol, optional
        The variable to solve for. If None, will use the first free symbol found.

    Returns:
    --------
    List[float]
        Sorted list of x-values where discontinuities occur within [x_min, x_max]

    Examples:
    ---------
    >>> find_discontinuities_in_range("1/(x-2)", 0, 5)
    [2.0]

    >>> find_discontinuities_in_range("1/((x-1)*(x+3))", -5, 5)
    [-3.0, 1.0]

    >>> find_discontinuities_in_range("1/((x-1)*(x+3))", 0, 5)
    [1.0]

    >>> find_discontinuities_in_range("(x+1)/(x**2 - 4)", -10, 10)
    [-2.0, 2.0]
    """
    """ # Convert string to sympy expression if needed
    if isinstance(expr, str):
        expr = sp.sympify(expr)

    # Access the individual terms using the .args attribute
    terms = expr.args
    print(f"Individual terms: {terms}") """

    # Determine the variable
    if var is None:
        free_vars = list(expr.free_symbols)
        if not free_vars:
            return []  # No variables, no discontinuities
        var = free_vars[0]
    elif isinstance(var, str):
        var = symbols(var)

    discontinuities = set()

    """ # 1. Find discontinuities from division by zero
    try:
        numer, denom = fraction(expr)
        if denom != 1:
            if numer.is_polynomial():
                numer = sp.factor(numer)

            if denom.is_polynomial():
                denom = sp.factor(denom)

            solution = sp.factor(numer/denom)
            numer, denom = solution.as_numer_denom()
            # Solve for where denominator equals zero
            zeros = solve(denom, var)
            for zero in zeros:
                if zero.is_real:
                    discontinuities.add(zero)
    except:
        pass """

    # 1. Find discontinuities from division by zero
    try:
        numer, denom = fraction(expr)
        if denom != 1:
            # Solve for where denominator equals zero
            zeros = solve_with_timeout(denom, var, sp.Interval(
                x_min, x_max, left_open=True, right_open=True))

            # zeros = solve(denom, var)
            for zero in zeros:
                if zero.is_real:
                    discontinuities.add(zero)
    except:
        pass

    # 2. Find discontinuities from logarithms
    for log_expr in expr.atoms(sp.log):
        arg = log_expr.args[0]
        # Log is undefined when argument <= 0
        critical_points = solve_with_timeout(arg, var, sp.Interval(
            x_min, x_max, left_open=True, right_open=True))
        # critical_points = solve(arg, var)
        for point in critical_points:
            if point.is_real:
                discontinuities.add(point)

    # 3. Find discontinuities from even roots (square root, 4th root, etc.)
    for pow_expr in expr.atoms(sp.Pow):
        base, exp = pow_expr.args
        # Check if exponent is a fraction with even denominator
        if exp.is_Rational and exp.q % 2 == 0 and exp.p > 0:
            # Even root is undefined for negative values
            try:
                # critical_points = sp.solve(base, var)

                critical_points = solve_with_timeout(base, var, sp.Interval(
                    x_min, x_max, left_open=True, right_open=True))

                if isinstance(critical_points, sp.FiniteSet):
                    for point in critical_points:
                        if point.is_real:
                            discontinuities.add(point)
            except:
                pass

    # 4. Filter discontinuities within the specified range [x_min, x_max]
    # filtered_discontinuities = []
    # for disc in discontinuities:
    #     try:
    #         #val = float(disc.evalf())
    #         # Check if within range (inclusive)
    #         #if x_min <= val <= x_max and abs(val) != float('inf'):
    #             #filtered_discontinuities.append(val)
    #     except:
    #         continue

    return discontinuities


def analyze_discontinuity_type(
    expr: Union[str, sp.Expr],
    x_val: float,
    var: Union[str, sp.Symbol] = None
) -> str:
    """
    Analyze the type of discontinuity at a specific point.

    Parameters:
    -----------
    expr : str or sympy expression
        The algebraic expression
    x_val : float
        The x-value to analyze
    var : str or sympy Symbol, optional
        The variable

    Returns:
    --------
    str
        Type of discontinuity: 'removable', 'infinite', 'jump', or 'unknown'
    """
    if isinstance(expr, str):
        expr = sp.sympify(expr)

    if var is None:
        free_vars = list(expr.free_symbols)
        if not free_vars:
            return 'unknown'
        var = free_vars[0]
    elif isinstance(var, str):
        var = symbols(var)

    # x_sym = symbols('x')
    # expr = expr.subs(var, x_sym)

    # Calculate left and right limits
    try:
        if x_val.is_real == False:
            return 'unknown'
        left_limit = sp.limit(expr, var, x_val, '-')
        right_limit = sp.limit(expr, var, x_val, '+')

        # if (expr.has(TrigonometricFunction)):
        #     # Infinite discontinuity
        #     if left_limit > 2e+14 or left_limit < -2e+14 or right_limit > 2e+14 or right_limit < -2e+14:
        #         return 'infinite'
        # else:
        # Infinite discontinuity
        if left_limit == oo or left_limit == -oo or right_limit == oo or right_limit == -oo:
            return 'infinite'

        # Removable discontinuity (limits equal but function undefined)
        if left_limit == right_limit and left_limit.is_finite:
            return ['removable', left_limit]

        # Jump discontinuity
        if left_limit != right_limit and left_limit.is_finite and right_limit.is_finite:
            return 'jump'

    except:
        pass

    return 'unknown'

# Unique list of elements from a given list


def unique_elements(input_list):
    unique_list = []
    for element in input_list:
        if element not in unique_list:
            unique_list.append(element)
    return unique_list


def pre_order_traversal(expression, detailed_results, x_min, x_max, var, level=0):
    for subtree in preorder_traversal(expression):
        if subtree.func.__name__ == "Mul" or subtree.func.__name__ == "Pow" or subtree.func.__name__ == "log":
            discontinuities = find_discontinuities_in_range(
                subtree, x_min, x_max, var)
            for disc in discontinuities:
                disc_type = analyze_discontinuity_type(full_expr, disc, var)
                if disc_type == 'unknown':
                    continue
                if isinstance(disc_type, list):
                    detailed_results.append([disc, disc_type[0], disc_type[1]])
                else:
                    detailed_results.append([disc, disc_type])


def find_discontinuities_detailed(
    expr: Union[str, sp.Expr],
    x_min: float,
    x_max: float,
    var: Union[str, sp.Symbol],
    period
) -> List[Tuple[float, str]]:
    """
    Find all discontinuities in an algebraic expression within a given range
    and classify their types.

    Parameters:
    -----------
    expr : str or sympy expression
        The algebraic expression to analyze
    x_min : float
        Lower bound of the range (inclusive)
    x_max : float
        Upper bound of the range (inclusive)
    var : str or sympy Symbol, optional
        The variable to solve for

    Returns:
    --------
    List[Tuple[float, str]]
        List of tuples (x_value, discontinuity_type) for each discontinuity

    Examples:
    ---------
    >>> find_discontinuities_detailed("1/(x-2)", 0, 5)
    [(2.0, 'infinite')]
    """

    # Convert string to sympy expression if needed
    if isinstance(var, str):
        var = symbols(var)

    # Convert string to sympy expression if needed
    if isinstance(expr, str):
        # sympify with evaluate=False to avoid side effects of sympy does not work
        # with trig functions that are overwrite.

        with sp.evaluate(False):
            expr = sp.sympify(expr)

    expr = trig_substitutions(expr)

    global full_expr
    full_expr = expr

    detailed_results = find_discontinuities(expr, var, x_min, x_max, period)

    # Remove duplicate discontinuities
    detailed_results = unique_elements(detailed_results)

    return detailed_results


app = Flask(__name__)


@app.route("/")
@app.route("/index")
def index():
    return render_template("index.html")


@app.route("/csolve", methods=['POST'])
def csolve():
    data = request.get_json()
    exp = data["exp"]
    var = data["var"]
    result = solve_for(exp, var)
    # print("solve")
    return jsonify({"result": result})


@app.route("/mode", methods=['POST'])
def mode():
    global mode_deg_rad
    data = request.get_json()
    m = data["mode"]
    # mode_deg_rad = m
    set_mode(m)
    return jsonify({"mode": m})


@app.route("/discontinuity", methods=['POST'])
def discontinuity():
    data = request.get_json()
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    if _exp == None:
        return []

    # while _exp.find("^") != -1:
    #     _exp = _exp.replace("^", "**")

    _exp = jsExpToPyExp(_exp)

    # disc = find_discontinuities(_exp, _var, lower, upper)
    # print(disc)

    # discont = discontinuities(_exp, lower, upper, _var)

    period = periodicity(sp.sympify(_exp), sp.Symbol(_var))
    if period != None:
        if get_mode() == "deg":
            period = period * 180 / sp.pi
        period = float(period)

    discont = find_discontinuities_detailed(
        _exp,
        lower, upper,
        _var,
        period
    )
    # discont = [0, "unknown2", 1]
    tps = []  # turning_points(_exp, _var, lower, upper)
    # print(discont)

    # for disc in discont:
    #     disc[0] = float(disc[0])
    #     if disc[2] != None:
    #         disc[2] = float(disc[2])

    # for disc in discont:
    #     disc[0] = float(disc[0])
    #     if len(disc) >= 3 and disc[2] != None:
    #         disc[2] = float(disc[2])

    # period = periodicity(sp.sympify(_exp), sp.Symbol(_var))
    # if period != None:
    #     if get_mode() == "deg":
    #         period = period * 180 / sp.pi
    #     period = float(period)

    return jsonify({"discontinuities": discont, "turningPoints": tps, "period": period})

# turning_points(expression, dependent_variable, lower_limit, upper_limit):


@app.route("/turningPoints", methods=['POST'])
def turningPoints():
    data = request.get_json()
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    _exp = jsExpToPyExp(_exp)

    turn_points = turning_points(_exp, lower, upper, _var)
    return jsonify({
        "turning_points": turn_points
    })


@app.route("/points", methods=['POST'])
def points():
    data = request.get_json()
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    # while _exp.find("^") != -1:
    #     _exp = _exp.replace("^", "**")

    _exp = jsExpToPyExp(_exp)

    # discont = discontinuities(_exp, lower, upper, _var)
    discont = find_discontinuities_in_range(
        _exp,
        lower, upper,
        _var
    )
    inflectn_points = inflection_points(_exp, lower, upper, _var)
    turn_points = turning_points(_exp, lower, upper, _var)
    return jsonify({
        "discontinuities": discont,
        "inflection_points": inflectn_points,
        "turning_points": turn_points
    })


if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=8080, threads=100)
