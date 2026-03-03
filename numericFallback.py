# %%fmt: off
import base64
import math
# import tempfile

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
from custom import custom
from sympy import lambdify
from sympy.functions.elementary.trigonometric import TrigonometricFunction
# from scipy.optimize import fsolve, root, root_scalar
import numpy as np
from sympy import symbols, solve, plot_implicit
import sympy as sp
# from domain_finder import closer_boundary
from degree_radian import sin_mode, cos_mode, tan_mode, cot_mode, sec_mode, csc_mode, asin_mode, acos_mode, atan_mode, acot_mode, asec_mode, acsc_mode, trig_substitutions
from my_misc import has_infinite_discontinuity_in_xrange, sanitize_contour_segments
import gc

from contourpy import contour_generator


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


# np.sqrt = custom_sqrt


# Define variables
x, y = symbols('x y')


def estimate_y_bounds2(equation, x_min, x_max, num_x=400, y_min=None, y_max=None, y_samples=400, match_tol=None, f_tol=1e-15):
    if (equation.has(TrigonometricFunction)) or ("_mode" in str(equation)):
        return (-300, 300)

    # Estimate y-range from symbolic solutions if possible
    x_vals = np.linspace(x_min, x_max, num_x)
    if y_min is None or y_max is None:
        try:
            y_sols = solve(equation, y)
        except Exception:
            y_sols = []

        if y_sols:
            y_vals_est = []
            for ys in y_sols:
                try:
                    y_fun = lambdify(x, ys, modules=[custom, 'numpy'])
                    y_eval = y_fun(x_vals)
                    y_eval = np.asarray(y_eval)
                    valid = ~np.isnan(y_eval) & np.isreal(y_eval)
                    if np.any(valid):
                        y_vals_est.append(np.real(y_eval[valid]))
                except Exception:
                    pass
            if y_vals_est:
                all_est = np.hstack(y_vals_est)
                est_min, est_max = float(
                    np.min(all_est)), float(np.max(all_est))
                padding = max(1.0, 0.1 * (est_max - est_min))
                if y_min is None:
                    y_min = est_min - padding
                if y_max is None:
                    y_max = est_max + padding

    # Fallback heuristic if still not set
    if y_min is None or y_max is None:
        y_guess = max(1.0, abs(x_min), abs(x_max)) * 10.0
        y_min = -y_guess if y_min is None else y_min
        y_max = y_guess if y_max is None else y_max

    return (y_min, y_max)


def generate_implicit_plot_points(expr, x_min=-10.0, x_max=10.0, has_discontinuity=False, y_min=-10.0, y_max=10.0):

    # if "sqrt" in str(expr):
    #     expr = str(expr).replace("sqrt", "custom_sqrt")
    #     expr = sp.sympify(expr)

    (_y_min, _y_max) = estimate_y_bounds2(expr, x_min, x_max)
    # (_y_min, _y_max) = (x_min, x_max)

    y_min = min(y_min, _y_min)
    y_max = max(y_max, _y_max)

    num_points = 800
    _x = np.linspace(x_min, x_max, num_points)
    _y = np.linspace(y_min, y_max, num_points)

    X, Y = np.meshgrid(_x, _y)
    # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
    f = lambdify((x, y), expr, modules=["numpy", custom])
    # z = z.astype(np.float32)
    z = f(X, Y)
    # z = z.astype(np.float32)
    z_val = 0.03*y_max
    # Convert the expression to a string
    expr_str = str(expr)
    if expr.has(TrigonometricFunction):
        z_val = np.maximum(z_val, 10)
    if ("_mode" in expr_str):
        z_val = np.maximum(z_val, 13)
    else:
        z_val = np.maximum(z_val, 15)
    # if abs(z_val) >= 1e100:
    # z_val = 10
    # z_val = np.percentile(z[~np.isnan(z)], 99)
    # z_val = np.nanmean(z)+3*np.nanstd(z)
    z_val = np.minimum(np.nanmax(z)*0.10, 45)
    # print(z_val)
    # z[np.abs(z) > z_val] = np.nan
    try:
        # Z = np.round(z, 2)  # Adjust precision as necessary
        # Replace inf with nan to avoid issues in contouring
        z[np.isinf(z)] = np.nan
        # CS = plt.contour(X, Y, np.ma.masked_where(z>z_val,
        #     z), levels=[0])
        # CS = plt.contour(X, Y, np.ma.masked_invalid(
        #     z), levels=[0], colors='blue', alpha=0)

        z_masked = np.ma.masked_where(z > z_val/4, z)

        cont_gen = contour_generator(
            X, Y, z=z_masked, name="serial")
        # cont_gen = contour_generator(X, Y, z, quad_as_tri=True, name="serial")
        del z
        del z_masked
        del X
        del Y
        del _x
        del _y
        gc.collect()  # Force garbage collection

        # lines(level) returns a list of branches (each is an (N, 2) array of coordinates)
        lines = cont_gen.lines(0)

        has_discontinuity = has_infinite_discontinuity_in_xrange(
            expr, x_min, x_max)

        all_points = []
        # for level_segments in CS.allsegs:
        # for level_segments in lines:
        for segment in lines:
            # segment is a NumPy array of shape (n_points, 2), where each row is [x, y]
            # all_points.append(segment)

            try:
                # Remove NaN values
                valid_mask = ~np.isnan(segment).any(axis=1)
                if not np.any(valid_mask):
                    segment = None  # Clear reference to segment to free memory
                    continue
                segment = segment[valid_mask]
            except Exception:
                # segment = None  # Clear reference to segment to free memory
                # continue
                pass

            # all_points.append(segment.astype(np.float32))
            segment = sanitize_contour_segments(
                expr, segment, x_min, x_max, has_discontinuity)
            if segment is None:
                continue
            all_points.append(base64.b64encode(
                segment.tobytes()).decode('utf-8'))
            # all_points.append(base64.b64encode(
            #     segment.astype(np.float32).tobytes()).decode('utf-8'))
            del segment

        del lines

        gc.collect()  # Force garbage collection
        return all_points, has_discontinuity

    except Exception as e:
        print(f"Error generating implicit plot points: {e}")
        return []


def get_sanitized_branches(expr, x_min, x_max, has_discontinuity, allsegs):
    for level_segments in allsegs:
        for branch in level_segments:
            branch = sanitize_contour_segments(
                expr, branch, x_min, x_max, has_discontinuity)
            yield branch.astype(np.float32)
            del branch
        del level_segments
        gc.collect()


# def generate_implicit_plot_points3(expr, x_min=-10.0, x_max=10.0, has_discontinuity=False, y_min=-10.0, y_max=10.0):

#     (_y_min, _y_max) = estimate_y_bounds2(expr, x_min, x_max)
#     # (_y_min, _y_max) = (x_min, x_max)

#     y_min = min(y_min, _y_min)
#     y_max = max(y_max, _y_max)

#     num_points = 1600
#     _x = np.linspace(x_min, x_max,np.round(num_points*0.88).astype(int))
#     _y = np.linspace(y_min, y_max, num_points)

#     X, Y = np.meshgrid(_x, _y)
#     # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
#     f = sp.lambdify((x, y), expr, modules=[custom, 'numpy'])
#     z = f(X, Y)
#     z_val = 0.03*y_max
#     # Convert the expression to a string
#     expr_str = str(expr)
#     if expr.has(TrigonometricFunction):
#         z_val = np.maximum(z_val, 10)
#     if ("_mode" in expr_str):
#         z_val = np.maximum(z_val, 13)
#     else:
#         z_val = np.maximum(z_val, 15)
#     # if abs(z_val) >= 1e100:
#     # z_val = 10
#     # z_val = np.percentile(z[~np.isnan(z)], 99)
#     # z_val = np.nanmean(z)+3*np.nanstd(z)
#     z_val = np.minimum(np.nanmax(z)*0.10, 45)
#     # print(z_val)
#     # z[np.abs(z) > z_val] = np.nan
#     try:
#         # Z = np.round(z, 2)  # Adjust precision as necessary
#         # Replace inf with nan to avoid issues in contouring
#         z[np.isinf(z)] = np.nan
#         cs = plt.contour(X, Y, np.ma.masked_where(z>z_val,
#             z), levels=[0])

#         # 3. Write to temporary file
#         line = 0
#         has_discontinuity = has_infinite_discontinuity_in_xrange(
#             expr, x_min, x_max)
#         with tempfile.NamedTemporaryFile( mode='w+',  delete=False, suffix=".txt") as temp_file:
#             for sanitized_branch in get_sanitized_branches(expr, x_min, x_max, has_discontinuity, cs.allsegs):
#                 if sanitized_branch is None:
#                     continue
#                 # Convert segment to bytes and base64
#                 encoded_branch = base64.b64encode(sanitized_branch.tobytes()).decode('utf-8')
#                 temp_file.write(encoded_branch + "\n")
#                 line += 1
#             temp_file.flush() # Ensure all data is written to disk
#             print(f"Branches written to: {temp_file.name} ({temp_file.tell()} bytes) ({line} lines)")
#             temp_file.close()

#         plt.clf()
#         plt.cla()
#         # gc.collect()
#         # plt.close()
#         # Mandatory cleanup

#         plt.close('all')
#         del cs
#         gc.collect()
#         return temp_file.name, has_discontinuity
#     except Exception as e:
#         print(f"Error generating implicit plot points: {e}")
#         return []

# def generate_implicit_plot_points4(expr, x_min=-10.0, x_max=10.0, has_discontinuity=False, y_min=-10.0, y_max=10.0):

#     (_y_min, _y_max) = estimate_y_bounds2(expr, x_min, x_max)
#     # (_y_min, _y_max) = (x_min, x_max)

#     y_min = min(y_min, _y_min)
#     y_max = max(y_max, _y_max)

#     num_points = 1600
#     _x = np.linspace(x_min, x_max,np.round(num_points*0.88).astype(int))
#     _y = np.linspace(y_min, y_max, num_points)

#     X, Y = np.meshgrid(_x, _y)
#     # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
#     f = sp.lambdify((x, y), expr, modules=[custom, 'numpy'])
#     z = f(X, Y)
#     z_val = 0.03*y_max
#     # Convert the expression to a string
#     expr_str = str(expr)
#     if expr.has(TrigonometricFunction):
#         z_val = np.maximum(z_val, 10)
#     if ("_mode" in expr_str):
#         z_val = np.maximum(z_val, 13)
#     else:
#         z_val = np.maximum(z_val, 15)
#     # if abs(z_val) >= 1e100:
#     # z_val = 10
#     # z_val = np.percentile(z[~np.isnan(z)], 99)
#     # z_val = np.nanmean(z)+3*np.nanstd(z)
#     z_val = np.minimum(np.nanmax(z)*0.10, 45)
#     # print(z_val)
#     # z[np.abs(z) > z_val] = np.nan
#     try:
#         # Z = np.round(z, 2)  # Adjust precision as necessary
#         # Replace inf with nan to avoid issues in contouring
#         z[np.isinf(z)] = np.nan
#         cs = plt.contour(X, Y, np.ma.masked_where(z>z_val,
#             z), levels=[0])

#         # 3. Write to temporary file
#         line = 0
#         has_discontinuity = has_infinite_discontinuity_in_xrange(
#             expr, x_min, x_max)
#         branches = []
#         for sanitized_branch in get_sanitized_branches(expr, x_min, x_max, has_discontinuity, cs.allsegs):
#             if sanitized_branch is None:
#                 continue
#             # Convert segment to bytes and base64
#             branches.append( base64.b64encode(sanitized_branch.tobytes()).decode('utf-8') )
#             line += 1

#         plt.clf()
#         plt.cla()
#         # gc.collect()
#         # plt.close()
#         # Mandatory cleanup

#         plt.close('all')
#         del cs
#         gc.collect()
#         return branches, has_discontinuity
#     except Exception as e:
#         print(f"Error generating implicit plot points: {e}")
#         return []
