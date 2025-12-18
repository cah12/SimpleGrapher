import sympy as sp


def find_extrema(expression, dependent_variable, lower_limit, upper_limit):
    # Create a SymPy symbol for the dependent variable
    var = sp.Symbol(dependent_variable)

    # Parse the expression using SymPy
    parsed_expression = sp.parse_expr(expression)

    # 1. Calculate the first derivative
    f_prime = sp.diff(parsed_expression, var)
    # print(f"First derivative: {f_prime}")

    if isinstance(f_prime, sp.Number):
        # print("The function is constant; no extrema exist.")
        return {
            'max': [],
            'min': []
        }

    # 2. Solve for critical points (where f_prime = 0)
    # solveset is generally preferred over solve
    critical_points = sp.solveset(f_prime, var, sp.Interval(
        lower_limit, upper_limit, left_open=True, right_open=True))
    # print(f"Critical points: {critical_points}")

    # 3. Calculate the second derivative
    f_double_prime = sp.diff(f_prime, var)
    # print(f"Second derivative: {f_double_prime}")

    # 4. Use the second derivative test
    max_points = []
    min_points = []
    num = 0
    for point in critical_points:
        if num > 200:
            # print("Too many critical points to evaluate.")
            break
        if point < lower_limit or point > upper_limit:
            num += 1
            continue
        # Check the sign of the second derivative at the critical point
        if f_double_prime.subs(var, point) < 0:
            max_x = point
            max_y = parsed_expression.subs(var, max_x)
            max_points.append([float(max_x), float(max_y)])
            # print(f"\nMaximum point (x, y): ({max_x}, {max_y})")
        elif f_double_prime.subs(var, point) > 0:
            min_x = point
            min_y = parsed_expression.subs(var, min_x)
            min_points.append([float(min_x), float(min_y)])
            # print(f"\nMinimum point (x, y): ({min_x}, {min_y})")
        else:
            # print(f"\nTest inconclusive for point {point}")
            pass

        num += 1

    # Return the maximum and minimum points as a dictionary
    return {
        'max': sorted(max_points, key=lambda x: x[0]),  # max_points,
        'min': sorted(min_points, key=lambda x: x[0])   # min_points
    }


expression = 'sin(4/x)'
dependent_variable = 'x'
lower_limit = -10
upper_limit = 10

result = find_extrema(expression, dependent_variable, lower_limit, upper_limit)
print(len(result['max']), len(result['min']))
# Output: {'max': 3, 'min': 0}
