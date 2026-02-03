import numpy as np
from scipy.optimize import brentq, newton
import warnings
import sympy as sp
from sympy import lambdify


def is_x_in_domain_numerical(f_xy, x_val, y_range=(-100, 100), tol=1e-15):
    """
    Checks if a given x_val is in the domain by trying to find a valid real y.

    Args:
        f_xy (function): The function f(x, y) = 0, which should return the value of the expression.
        x_val (float): The specific x value to check.
        y_range (tuple): The range (min, max) to search for a valid y.
        tol (float): Tolerance for the root finder.

    Returns:
        bool: True if a valid real y is found for the given x, False otherwise.
    """
    def f_y(y):
        # Create a partial function for y given a fixed x
        return f_xy(x_val, y)

    try:
        # We look for a root (a value of y where f(x_val, y) == 0) within the y_range.
        # This requires f(y_min) and f(y_max) to have opposite signs, which is a limitation.
        # A more robust check is to iterate/search for potential roots.

        # A common numerical way to check existence is to look for sign changes in a range
        # However, a simpler, more direct check for *validity* (e.g. no NaNs) is possible.

        # Let's check for basic numerical validity without assuming sign change first.
        # The primary conditions for a valid domain in numerical contexts are often
        # related to avoiding undefined operations (like sqrt of negative, division by zero)

        # For a general function, we can sample the y range and use a root finder if sign changes occur.

        # A simpler check: Does evaluating the function produce a real number?
        # This doesn't guarantee f(x,y) = 0 can be met, but it checks if f is defined for any y.
        # The user's goal is to find x for which a *solution* y exists.

        # Numerical root finding for y given x
        # Sample points in y_range to find an initial bracket
        y_samples = np.linspace(y_range[0], y_range[1], 100)
        signs = np.sign(f_y(y_samples))
        for i in range(len(signs) - 1):
            if signs[i] * signs[i+1] < 0:
                # Sign change found, a root exists in this small interval
                # Use a root finder for more precision
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    brentq(f_y, y_samples[i], y_samples[i+1], xtol=tol)
                return True  # A valid y was found

        # Also check near zero if the range includes it, for cases where f(y) might touch zero without crossing
        if any(np.abs(f_y(y_samples)) < tol):
            return True

        return False  # No sign change or near-zero value found

    except (ValueError, ZeroDivisionError, TypeError, FloatingPointError):
        # Catches common numerical errors (e.g., sqrt(-1), division by zero, invalid input types)
        return False

# Example Usage:
# Equation: x**2 + y**2 - 4 = 0 (A circle with radius 2. Domain for x is [-2, 2])


# def circle_func(x, y):
#     # return x**2 + y**2 - 4
#     return y**6+8*y-x
# # x, y = sp.symbols('x y')
# # exp = x**2 + y**2 - 4
# # circle_func = lambdify([x, y], exp)


# # Test cases
# x_test1 = 1.0  # Should be in domain
# x_test2 = 3.0  # Should be outside domain
# x_test3 = -2.0  # Should be in domain

# current_boundary = (-6.784, -0.9276)
# next_point = (-6.683, -0.9034)
def closer_boundary(fn, current_boundary, next_point):
    y_step = np.abs(current_boundary[1] - next_point[1])
    x_step = np.abs(current_boundary[0] - next_point[0])

    # y_step = np.max([y_step, 0.01])
    # x_step = np.max([x_step, 0.01])

    y_range = (current_boundary[0]-6*y_step, current_boundary[1]+6*y_step)
    x_range = (current_boundary[0]-6*x_step, current_boundary[1]+6*x_step)

    x_samples = np.linspace(x_range[0], x_range[1], 500)
    y_samples = np.linspace(y_range[0], y_range[1], 1000)

    for x in x_samples:
        # print(x, is_x_in_domain_numerical(circle_func, x, y_range))
        if is_x_in_domain_numerical(fn, x, y_range):
            y = sp.symbols('y')
            _fn = fn(x, y)

            for i in range(len(y_samples) - 1):
                if _fn.subs(y, y_samples[i]) * _fn.subs(y, y_samples[i+1]) < 0:
                    return (x.item(), y_samples[i].item())

    return None


# first point on left boundary or last point on right boundary
# current_boundary = (-6.784, -0.9276)
# next_point = (-6.683, -0.9034)

# print(closer_boundary(circle_func, current_boundary, next_point))

# print(
#     f"Is x = {x_test1} in the domain? {is_x_in_domain_numerical(circle_func, x_test1)}")
# print(
#     f"Is x = {x_test2} in the domain? {is_x_in_domain_numerical(circle_func, x_test2)}")
# print(
#     f"Is x = {x_test3} in the domain? {is_x_in_domain_numerical(circle_func, x_test3)}")
