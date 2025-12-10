from sympy import symbols, singularities, solveset, cos, sin, pi
from sympy.calculus.util import continuous_domain, Interval

x = symbols('x')

# Example 1: A rational function
# expr_rational = 1 / (x**2 - 1)
expr_rational = 1 / sin(x)
sings_rational = singularities(expr_rational, x)
print(f"Singularities of {expr_rational}: {sings_rational}")
# Output: Singularities of 1/(x**2 - 1): {-1, 1}

# # Example 2: A function with trigonometric singularities
# expr_trig = 1 / cos(x)
# # For trigonometric functions, you often need to specify a domain
# domain = Interval(0, 2*pi)
# # Note: The singularities function works well for specific domains.
# # For general cases without a domain, SymPy might return a ConditionSet.

# # A common approach for rational functions (like 1/cos(x)) is to solve the denominator:
# denominator = cos(x)
# sings_trig = solveset(denominator, x, domain)
# print(f"Singularities of {expr_trig} in [0, 2pi]: {sings_trig}")
# # Output: Singularities of 1/cos(x) in [0, 2pi]: {pi/2, 3*pi/2}
