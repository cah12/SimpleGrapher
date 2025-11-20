from ast import expr
from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

import sympy as sp
from sympy import symbols, solve, fraction, oo, S
from typing import List, Union, Tuple


# The Degree To Radian Code
# --------------------------#
####################################################
###########################################################
mode_deg_rad = "deg"


def deg2rad(d):
    return sp.N(d*sp.pi/180)


original_sin = sp.sin


def sin(d):
    if mode_deg_rad == "deg":
        return original_sin(deg2rad(d))
    return original_sin(d)


sp.sin = sin

original_cos = sp.cos


def cos(d):
    if mode_deg_rad == "deg":
        return original_cos(deg2rad(d))
    return original_cos(d)


sp.cos = cos

original_tan = sp.tan


def tan(d):
    if mode_deg_rad == "deg":
        return original_tan(deg2rad(d))
    return original_tan(d)


sp.tan = tan

original_sec = sp.sec


def sec(d):
    if mode_deg_rad == "deg":
        return 1/original_cos(deg2rad(d))
    return 1/original_cos(d)


sp.sec = sec

original_csc = sp.csc


def csc(d):
    if mode_deg_rad == "deg":
        return 1/original_sin(deg2rad(d))
    return 1/original_sin(d)


sp.csc = csc

original_cot = sp.cot


def cot(d):
    if mode_deg_rad == "deg":
        return 1/original_tan(deg2rad(d))
    return 1/original_tan(d)


sp.cot = cot


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
        if arr[0] == c:
            return [pyExpToJsExp(arr[1])]
        if arr[1] == c:
            return [pyExpToJsExp(arr[0])]

    if len(arr) == 1:
        arr.append("0")

    result = []
    try:
        if arr[1] == "0":
            solutions = sp.solveset(sp.parse_expr(arr[0]), v)
        else:
            s = sp.solveset(
                sp.Eq(sp.parse_expr(arr[0]), sp.parse_expr(arr[1])), v)
            solutions = list(s)

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


def turning_points(expr, lower, upper, var):
    try:
        d = sp.diff(sp.parse_expr(expr), sp.Symbol(var))
        dis = sp.solveset(d, var, sp.Interval(
            lower, upper, left_open=True, right_open=True))
        return list(map(float, dis))
    except:
        return []


def inflection_points(expr, lower, upper, var):
    try:
        x = sp.Symbol(var)
        d = sp.diff(sp.parse_expr(expr), x, x)
        dis = sp.solveset(d, var, sp.Interval(lower, upper))
        return list(map(float, dis))
    except:
        return []


""" def discontinuities(exp_, lower, upper, _var):  
    discount = []  
    v = sp.Symbol(_var)

    _exp=sp.parse_expr(exp_) 
    
    num, denom = _exp.as_numer_denom()  

    if denom == 1:
        return []  
    
    if num.is_polynomial():
        num = sp.factor(num)  

    if denom.is_polynomial():
        denom = sp.factor(denom)             

    
    solution = sp.factor(num/denom)
    num, denom = solution.as_numer_denom()
    ds = sp.solveset(denom, v, sp.Interval(lower, upper, left_open=True, right_open=True))   
    

    if type(ds)==sp.FiniteSet:
        for sol in list(ds):
            try:
                v = float(sol.evalf())
                discount.append(v)
            except:
                pass
  
    
    if len(discount) == 0:
        try:
            d= list(sp.singularities(_exp, v))
            for sol in d:
                try:
                    v = float(sol.evalf())
                    discount.append(v)
                except:
                    pass

            if len(discount):
                arr = []
                for s in discount:
                    if s < upper and s > lower:
                        arr.append(s)
                discount = arr
        except:
            pass

            

    
    discount.sort()   
    discount = list(map(float, discount))
    
    return discount    """


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
            zeros = sp.solveset(denom, var, sp.Interval(
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
        critical_points = sp.solveset(arg, var, sp.Interval(
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
            critical_points = sp.solveset(base, var, sp.Interval(
                x_min, x_max, left_open=True, right_open=True))
            # critical_points = solve(base, var)
            for point in critical_points:
                if point.is_real:
                    discontinuities.add(point)

    # 4. Filter discontinuities within the specified range [x_min, x_max]
    filtered_discontinuities = []
    for disc in discontinuities:
        try:
            val = float(disc.evalf())
            # Check if within range (inclusive)
            if x_min <= val <= x_max and abs(val) != float('inf'):
                filtered_discontinuities.append(val)
        except:
            continue

    return sorted(filtered_discontinuities)


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

    x_sym = symbols('x')
    expr = expr.subs(var, x_sym)

    # Calculate left and right limits
    try:
        left_limit = sp.limit(expr, x_sym, x_val, '-')
        right_limit = sp.limit(expr, x_sym, x_val, '+')

        # Infinite discontinuity
        if left_limit == oo or left_limit == -oo or right_limit == oo or right_limit == -oo:
            return 'infinite'

        # Removable discontinuity (limits equal but function undefined)
        if left_limit == right_limit and left_limit.is_finite:
            return ['removable', float(left_limit.evalf())]

        # Jump discontinuity
        if left_limit != right_limit and left_limit.is_finite and right_limit.is_finite:
            return 'jump'

    except:
        pass

    return 'unknown'


def find_discontinuities_detailed(
    expr: Union[str, sp.Expr],
    x_min: float,
    x_max: float,
    var: Union[str, sp.Symbol] = None
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
    if isinstance(expr, str):
        expr = sp.sympify(expr)

    # Get the terms
    terms = sp.Add.make_args(expr)

    # print("Expression:", expr)
    # print("Terms:", terms)

    # discontinuities = find_discontinuities_in_range(expr, x_min, x_max, var)

    detailed_results = []
    for term in terms:
        discontinuities = find_discontinuities_in_range(
            term, x_min, x_max, var)
        for disc in discontinuities:
            disc_type = analyze_discontinuity_type(expr, disc, var)
            if isinstance(disc_type, list):
                detailed_results.append([disc, disc_type[0], disc_type[1]])
            else:
                detailed_results.append([disc, disc_type])

    # print("Discontinuities:", detailed_results)

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
    mode_deg_rad = m
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

    # print("discontinuity")
    # discont = discontinuities(_exp, lower, upper, _var)
    discont = find_discontinuities_detailed(
        _exp,
        lower, upper,
        _var
    )
    # print(discont)

    return jsonify({"discontinuities": discont})


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
