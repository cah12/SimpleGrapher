
# import math

import numpy as np
import sympy as sp


# def custom_sqrt_func(x):
#     """Calculates the square root of a single value.
#     Handles negative numbers by returning a complex number."""
#     if x >= 0:
#         return math.sqrt(x)
#     else:
#         return math.sqrt(abs(x)) * -1


# # custom_sqrt_ufunc is now a universal function
# custom_sqrt_ufunc = np.frompyfunc(custom_sqrt_func, 1, 1)

# np.sqrt = custom_sqrt_ufunc

# def custom_sqrt(x):
#     """
#     Calculates sqrt(abs(x)) for non-negative x, and -sqrt(abs(x)) for negative x.
#     """
#     # Convert input to a NumPy array for element-wise operations if it isn't already one
#     x = np.array(x)

#     # Condition: elements less than 0
#     is_negative = x < 0

#     # Calculate the desired values for both cases
#     # For all elements, calculate -sqrt(abs(x))
#     negative_result = -np.sqrt(np.abs(x))
#     # For positive elements, calculate sqrt(x) which is also sqrt(abs(x))
#     positive_result = np.sqrt(x)

#     # Use np.where to select from the two results based on the condition
#     # If is_negative is True, use the value from negative_result
#     # Otherwise, use the value from positive_result
#     result = np.where(is_negative, negative_result, positive_result)
#     return result


def custom_sin_mode(arg):
    return np.sin(np.deg2rad(arg))


def custom_cos_mode(arg):
    return np.cos(np.deg2rad(arg))


def custom_tan_mode(arg):
    return np.tan(np.deg2rad(arg))


def custom_cot_mode(arg):
    return 1 / np.tan(np.deg2rad(arg))


def custom_sec_mode(arg):
    return 1 / np.cos(np.deg2rad(arg))


def custom_csc_mode(arg):
    return 1 / np.sin(np.deg2rad(arg))


def custom_asin_mode(arg):
    return np.rad2deg(np.arcsin(arg))


def custom_acos_mode(arg):
    return np.rad2deg(np.arccos(arg))


def custom_atan_mode(arg):
    return np.rad2deg(np.arctan(arg))


def custom_acot_mode(arg):
    return np.rad2deg(np.arccot(arg))


def custom_asec_mode(arg):
    return np.rad2deg(np.arccos(1/arg))


def custom_acsc_mode(arg):
    return np.rad2deg(np.arcsin(1/arg))


custom = {
    "sin_mode": custom_sin_mode,
    "cos_mode": custom_cos_mode,
    "tan_mode": custom_tan_mode,
    "cot_mode": custom_cot_mode,
    "sec_mode": custom_sec_mode,
    "csc_mode": custom_csc_mode,
    "asin_mode": custom_asin_mode,
    "acos_mode": custom_acos_mode,
    "atan_mode": custom_atan_mode,
    "acot_mode": custom_acot_mode,
    "asec_mode": custom_asec_mode,
    "acsc_mode": custom_acsc_mode
}
