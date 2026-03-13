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
import re

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
    # if (equation.has(TrigonometricFunction)) or ("_mode" in str(equation)):
    #     return (-300, 300)

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


# def estimate_y_bounds2(equation, x_min, x_max, num_x=400, y_min=None, y_max=None, y_samples=400, match_tol=None, f_tol=1e-15):
#     if (equation.has(TrigonometricFunction)) or ("_mode" in str(equation)):
#         return (-300, 300)

#     # Estimate y-range from symbolic solutions if possible
#     x_vals = np.linspace(x_min, x_max, num_x)
#     if y_min is None or y_max is None:
#         try:
#             y_sols = solve(equation, y)
#         except Exception:
#             y_sols = []

#         if y_sols:
#             y_vals_est = []
#             for ys in y_sols:
#                 try:
#                     y_fun = lambdify(x, ys, modules=[custom, 'numpy'])
#                     y_eval = y_fun(x_vals)
#                     y_eval = np.asarray(y_eval)
#                     valid = ~np.isnan(y_eval) & np.isreal(y_eval)
#                     if np.any(valid):
#                         y_vals_est.append(np.real(y_eval[valid]))
#                 except Exception:
#                     pass
#             if y_vals_est:
#                 all_est = np.hstack(y_vals_est)
#                 est_min, est_max = float(
#                     np.min(all_est)), float(np.max(all_est))
#                 padding = max(1.0, 0.1 * (est_max - est_min))
#                 if y_min is None:
#                     y_min = est_min - padding
#                 if y_max is None:
#                     y_max = est_max + padding

#     # Fallback heuristic if still not set
#     if y_min is None or y_max is None:
#         y_guess = max(1.0, abs(x_min), abs(x_max)) * 10.0
#         y_min = -y_guess if y_min is None else y_min
#         y_max = y_guess if y_max is None else y_max

#     return (y_min, y_max)


def grid_x_y_z_val(expr, x_min, x_max, y_min, y_max):
    z_val = 8
    has_discontinuity = has_infinite_discontinuity_in_xrange(
        expr, x_min, x_max)

    num_x = 200
    density = 1000000

    # if has_discontinuity or "/sin" in str(expr) or "/cos" in str(expr):
    #     num_x = 1000
    # density = 1000000

    num_y_max = 1000
    num_y = num_y_max

    excpt = False

    # Pattern for common trig functions
    pattern = r"\b(sin|cos|tan|asin|acos|atan|sin_mode|cos_mode|tan_mode|asin_mode|acos_mode|atan_mode)\b"

    # Replace all matched trig functions with '1'
    str_expr = re.sub(pattern, "0*", str(expr))
    str_expr = sp.sympify(str_expr)

    d_x = 1

    if str_expr.is_polynomial(x):
        d_x = int(sp.degree(str_expr, x))

    d_y = 1

    if str_expr.is_polynomial(y):
        d_y = int(sp.degree(str_expr, y))

    try:
        # num_y = np.minimum(num_y_max, int(num_x*(y_max-y_min)/(x_max-x_min)))
        num_y = int(num_x*(y_max-y_min)/(x_max-x_min))
    except Exception:
        excpt = True
        pass

    if excpt == True:
        num_x = int(num_y_max*2)
        num_y = num_x
        # num_x = int(density/(num_y_max*2))

    while num_y*num_x < density/50:
        num_x *= 2
        num_y = int(num_x*(y_max-y_min)/(x_max-x_min))

    # if has_discontinuity or "/sin" in str(expr) or "/cos" in str(expr):
    if has_discontinuity:
        # if "/sin" in str(expr) or "/cos" in str(expr):
        # if expr has 1 y and degree of poly in x is less than 3
        # num_x = 1000
        # num_y = 1000
        z_val = 8

        # if expr has 1 y and degree of poly in x is greater than 2
        if expr.has(TrigonometricFunction):
            # num_x = 400
            # num_y = 400
            d_y = np.minimum(d_y, 6)
            if d_x == 1:
                z_val = 8*d_y
            elif d_y == 1:
                z_val = 8*d_x
            else:
                z_val = 8*d_x*d_y

            z_val = 10
        # if "_mode" in str(expr):
        #     z_val = 15

        # return num_x, num_y, z_val, has_discontinuity
    elif (expr.has(TrigonometricFunction) or "_mode" in str(expr)):
        # num_x = 200
        # num_y = 200
        z_val = 8
        # return num_x, num_y, z_val, has_discontinuity

    else:
        # x^3 = y, x^3-40 = y, x^3-40x^2 = y, x^3-40x+10 = y
        # num_x = 400
        # num_y = 400
        z_val = np.finfo(np.float64).max

    # d_y = 1
    # if expr.is_polynomial(y):
    #     d_y = sp.degree(expr, y)
    # if 5 >= d_y >= 3:
    #     num_y = np.maximum(num_y, 1000)
    # elif 30 >= d_y > 5:
    #     num_y = np.maximum(num_y, 2500)
    # elif d_y > 30:
    #     num_y = np.maximum(num_y, 3000)
    #     # num_x = np.maximum(num_x, 2000)

    # if num_y <= num_x and d_y > 10:
    #     num_y = int(2*num_x)

    # return num_x, num_y, z_val, has_discontinuity
    factor = 1
    default_value = 50
    if 7 >= d_x >= 3:
        factor = 2
    elif 15 >= d_x > 7:
        factor = 3
    elif 20 >= d_x > 15:
        factor = 4
    elif 30 >= d_x > 20:
        factor = 5
    elif 40 >= d_x > 30:
        factor = 6
    elif 50 >= d_x > 40:
        factor = 7
    elif d_x > 50:
        factor = 8

    factor_y = 1

    if 7 >= d_y >= 3:
        factor_y = 2
    elif 15 >= d_y > 7:
        factor_y = 3
    elif 20 >= d_y > 15:
        factor_y = 5
    elif 30 >= d_y > 20:
        factor_y = 6
    elif 40 >= d_y > 30:
        factor_y = 7
    elif 50 >= d_y > 40:
        factor_y = 8
    elif d_y > 50:
        factor_y = 9

    if has_discontinuity == False and expr.has(TrigonometricFunction):
        default_value *= 10
    if has_discontinuity:
        default_value *= 10

    if d_x > 3:
        z_val = np.finfo(np.float64).max
    return default_value*factor_y*2, default_value*factor_y*2, z_val, has_discontinuity


# class QuadTreeNode:
#     def __init__(self, value, x_start, x_end, y_start, y_end):
#         self.value = value
#         self.x_start = x_start
#         self.x_end = x_end
#         self.y_start = y_start
#         self.y_end = y_end
#         self.children = [None, None, None, None]

#     def split(self):
#         x_mid = (self.x_start + self.x_end) / 2
#         y_mid = (self.y_start + self.y_end) / 2
#         self.children[0] = QuadTreeNode(
#             self.value, self.x_start, x_mid, self.y_start, y_mid)
#         self.children[1] = QuadTreeNode(
#             self.value, x_mid, self.x_end, self.y_start, y_mid)
#         self.children[2] = QuadTreeNode(
#             self.value, self.x_start, x_mid, y_mid, self.y_end)
#         self.children[3] = QuadTreeNode(
#             self.value, x_mid, self.x_end, y_mid, self.y_end)

#     def update(self, value):
#         self.value = value


# def adaptive_patch(f, Z, patch_size, threshold=0.1):
#     root = QuadTreeNode(np.mean(Z), 0, 1, 0, 1)
#     root.split()
#     queue = [root]
#     Z_array = np.zeros((patch_size, patch_size))

#     while queue:
#         node = queue.pop(0)
#         if node.value is not None:
#             cell_vals = [
#                 node.value] + [child.value for child in node.children if child is not None]
#             if min(cell_vals) < threshold and max(cell_vals) > -threshold:
#                 x_start, x_end = node.x_start, node.x_end
#                 y_start, y_end = node.y_start, node.y_end
#                 x_sub = np.linspace(x_start, x_end, patch_size)
#                 y_sub = np.linspace(y_start, y_end, patch_size)
#                 sub_patch = f(x_sub, y_sub)
#                 node.update(np.mean(sub_patch))
#                 for i, child in enumerate(node.children):
#                     if child is None:
#                         child = QuadTreeNode(np.mean(Z[node.y_start:node.y_end, node.x_start:node.x_end]),
#                                              node.x_start + i/4 *
#                                              (node.x_end - node.x_start),
#                                              node.x_start +
#                                              (i+1)/4 * (node.x_end - node.x_start),
#                                              node.y_start + i/4 *
#                                              (node.y_end - node.y_start),
#                                              node.y_start + (i+1)/4 * (node.y_end - node.y_start))
#                         node.children[i] = child
#                     queue.append(child)
#                 Z[node.y_start:node.y_end, node.x_start:node.x_end] = sub_patch
#         else:
#             queue.insert(0, node)
#     return Z

def adaptive_patch(f, Z, patch_size, n_base_x, n_base_y, x_coarse, y_coarse, threshold=0.01):
    for i in range(0, n_base_y - 1):
        for j in range(0, n_base_x - 1):
            # Check if this cell boundary crosses 0
            cell_vals = [Z[i, j], Z[i+1, j], Z[i, j+1], Z[i+1, j+1]]
            if min(cell_vals) < threshold and max(cell_vals) > -threshold:
                # Calculate physical boundaries for this cell
                x_start, x_end = x_coarse[j], x_coarse[j+1]
                y_start, y_end = y_coarse[i], y_coarse[i+1]

                # 3. Create High-Resolution Subarray
                x_sub = np.linspace(x_start, x_end, patch_size)
                y_sub = np.linspace(y_start, y_end, patch_size)
                xx, yy = np.meshgrid(x_sub, y_sub)
                sub_patch = f(xx, yy)

                # 4. Patch the main array (using coarse indices to map)
                # This example demonstrates placing the patch.
                # Real adaptive mesh often uses a separate data structure.
                # To keep z as one array, we must map patch_size to the
                # pixel density of the main array.

                # Dummy placement for illustration (usually handled by quadtrees)
                # For simplicity, we just update a region of the coarse mesh
                # with the average of the patch here.
                # z[i:i+4, j:j+4] = patch
                if i < n_base_y - patch_size and j < n_base_x - patch_size:
                    Z[i:i+patch_size, j:j+patch_size] = sub_patch
    return Z


def generate_implicit_plot_points(expr, x_min=-10.0, x_max=10.0, autoScale=False, has_discontinuity=False, y_min=-10.0, y_max=10.0):

    # if "sqrt" in str(expr):
    # expr = str(expr).replace("sin", "np.sin")
    # expr = sp.sympify(expr)

    # width = x_max - x_min
    # x_min = x_min - 1 * width
    # x_max = x_max + 1 * width

    if autoScale:
        (y_min, y_max) = estimate_y_bounds2(expr, x_min, x_max)
        # (_y_min, _y_max) = (x_min, x_max)

        # y_min = min(y_min, _y_min)
        # y_max = max(y_max, _y_max)

        _max = np.max([np.abs(y_min), np.abs(y_max)])

        # if expr.has(TrigonometricFunction):
        # _max = np.minimum(_max, 1e16)

        y_min = -_max
        y_max = _max

    # if the power of y > 1
    # d_p = sp.degree(expr, gen=y

    # if the power of y == 1
    num_points = 800
    d = 4

    # d_p = sp.degree(expr, gen=y)
    # # if the power of y > 1
    # if d_p > 1:
    #     num_points = 2000
    #     d = 1

    num_x, num_y, z_val, has_discontinuity = grid_x_y_z_val(
        expr, x_min, x_max, y_min, y_max)

    # num_x = 500
    # num_y = 500

    _x = np.linspace(x_min, x_max, num_x)
    _y = np.linspace(y_min, y_max, num_y)
    # _x = np.linspace(x_min, x_max, num_points)
    # _y = np.linspace(y_min, y_max, num_points)

    X, Y = np.meshgrid(_x, _y)
    # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
    f = lambdify((x, y), expr, modules=["numpy", custom])
    # z = z.astype(np.float32)
    z = f(X, Y)

    # z = z.astype(np.float32)
    # z_val = 0.03*y_max
    # # Convert the expression to a string
    # expr_str = str(expr)
    # if expr.has(TrigonometricFunction):
    #     z_val = np.maximum(z_val, 10)
    # if ("_mode" in expr_str):
    #     z_val = np.maximum(z_val, 13)
    # else:
    #     z_val = np.maximum(z_val, 15)
    # # if abs(z_val) >= 1e100:
    # # z_val = 10
    # # z_val = np.percentile(z[~np.isnan(z)], 99)
    # # z_val = np.nanmean(z)+3*np.nanstd(z)
    # z_val = np.minimum(np.nanmax(z)*0.10, 45)
    # print(z_val)
    # z[np.abs(z) > z_val] = np.nan
    large_range_span = False
    try:
        # Z = np.round(z, 2)  # Adjust precision as necessary
        # Replace inf with nan to avoid issues in contouring
        z[np.isinf(z)] = np.nan
        # CS = plt.contour(X, Y, np.ma.masked_where(z>z_val,
        #     z), levels=[0])
        # CS = plt.contour(X, Y, np.ma.masked_invalid(
        #     z), levels=[0], colors='blue', alpha=0)

        # z = adaptive_patch(f, z, 10, num_x, num_y, _x, _y)
        # z = adaptive_patch(f, z, 4)

        z_masked = np.ma.masked_where(z < -z_val, z)

        # z[np.abs(z) < 1e-4] = 0.0

        cont_gen = contour_generator(
            X, Y, z=z_masked, name="serial")
        # cont_gen = contour_generator(X, Y, z, quad_as_tri=True, corner_as_point=True, name="serial")
        del z
        del z_masked
        del X
        del Y
        del _x
        del _y
        # gc.collect()  # Force garbage collection

        # lines(level) returns a list of branches (each is an (N, 2) array of coordinates)
        lines = cont_gen.lines(0)

        cont_gen = None
        gc.collect()

        # has_discontinuity = has_infinite_discontinuity_in_xrange(
        #     expr, x_min, x_max)

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
            if not large_range_span:
                max_y = np.max(segment[:, 1])
                min_y = np.min(segment[:, 1])
                if np.abs(max_y - min_y) > 1e16:
                    large_range_span = True

            all_points.append(base64.b64encode(
                segment.tobytes()).decode('utf-8'))
            # all_points.append(base64.b64encode(
            #     segment.astype(np.float32).tobytes()).decode('utf-8'))
            del segment

        del lines

        gc.collect()  # Force garbage collection
        return all_points, has_discontinuity, large_range_span

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
