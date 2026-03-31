# %%fmt: off
import base64
from custom import custom
from sympy import lambdify
from sympy.functions.elementary.trigonometric import TrigonometricFunction
import numpy as np
from sympy import symbols, solve
import sympy as sp
from my_misc import has_infinite_discontinuity_in_xrange, sanitize_contour_segments
import re
from contourpy import contour_generator

# Enable saving of uncollectable objects
# gc.set_debug(gc.DEBUG_SAVEALL)


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
                # if y_min is None:
                #     y_min = est_min - padding
                # if y_max is None:
                #     y_max = est_max + padding
                if y_min is None:
                    y_min = est_min - 0.1*abs(est_min)
                if y_max is None:
                    y_max = est_max + 0.1*abs(est_max)

    # Fallback heuristic if still not set
    if y_min is None or y_max is None:
        y_guess = max(1.0, abs(x_min), abs(x_max)) * 10.0
        y_min = -y_guess if y_min is None else y_min
        y_max = y_guess if y_max is None else y_max
    x_vals = None
    return (y_min, y_max)


def grid_x_y_z_val(expr, x_min, x_max, y_min, y_max):
    z_val = 8
    has_discontinuity = has_infinite_discontinuity_in_xrange(
        expr, x_min, x_max)

    num_x = 200
    density = 1000000

    # if has_discontinuity or "/sin" in str(expr) or "/cos" in str(expr):
    #     num_x = 1000
    # density = 1000000

    num_y_max = 500
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
            if d_x == 0:
                pass
            elif d_x == 1:
                z_val = 8.0*d_y
            elif d_y == 1:
                z_val = 8*(d_x+0.8)
            else:
                z_val = 8*(d_x+0.8)*d_y

            # z_val = 10
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
        if not (expr.has(TrigonometricFunction) or "_mode" in str(expr)) or d_x > 10:
            z_val = np.finfo(np.float64).max
        else:
            z_val = z_val**d_x
    return default_value*factor_y, default_value*factor_y, z_val, has_discontinuity


# def estimate_z_val111(expr, x_min, x_max, y_min, y_max):
#     """
#     Estimate a reasonable z_val for an implicit function given the domain (x_min, x_max).

#     This function computes the maximum absolute value of the function over the domain and
#     returns 10% of that value as the estimated z_val. If the maximum absolute value is
#     less than 1e-100, then the function returns 1e-100.

#     Args:
#         expr (sympy expression): The implicit function to estimate z_val for.
#         x_min (float): The minimum x-value of the domain.
#         x_max (float): The maximum x-value of the domain.

#     Returns:
#         float: The estimated z_val value.
#     """

#     # Compute the maximum absolute value of the function over the domain
#     x = sp.Symbol('x')
#     f = sp.lambdify((x, y), expr, modules=['numpy'])
#     x_val = np.linspace(x_min, x_max, 500)
#     y_val = np.linspace(y_min, y_max, 500)
#     max_abs_val = np.max(np.abs(f(x_val, y_val)))

#     # If the maximum absolute value is less than 1e-100, return 1e-100
#     if max_abs_val < 1e-100:
#         return 1e-100

#     # Otherwise, return 10% of the maximum absolute value as the estimated z_val
#     return max_abs_val * 0.08


def generate_implicit_plot_points(expr, x_min=-10.0, x_max=10.0, autoScale=False, has_discontinuity=False, y_min=-10.0, y_max=10.0):

    # if "sqrt" in str(expr):
    # expr = str(expr).replace("sin", "np.sin")
    # expr = sp.sympify(expr)

    # width = x_max - x_min
    # x_min = x_min - 1 * width
    # x_max = x_max + 1 * width

    # if autoScale:
    (y_min, y_max) = estimate_y_bounds2(expr, x_min, x_max)

    # (_y_min, _y_max) = (x_min, x_max)

    # y_min = min(y_min, _y_min)
    # y_max = max(y_max, _y_max)

    _max = np.max([np.abs(y_min), np.abs(y_max)])

    # if expr.has(TrigonometricFunction):
    # _max = np.minimum(_max, 1e16)

    y_min = -_max
    y_max = _max

    # y_min = -100
    # y_max = 100

    # print(y_min, y_max)

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

    num_x = 500
    num_y = 500

    _x = np.linspace(x_min, x_max, num_x)
    if y_max > 1e16:
        f = 0.095
        _y = np.linspace(y_min, -1e12, int(num_y*f), endpoint=False)
        _y = np.append(_y, np.linspace(-1e12, 1e12,
                       int(num_y*(1-2*f)), endpoint=False))
        _y = np.append(_y, np.linspace(1e12, y_max, int(num_y*f)))
    else:
        _y = np.linspace(y_min, y_max, num_y)

    X, Y = np.meshgrid(_x, _y)

    f = lambdify((x, y), expr, modules=["numpy", custom])
    # z = z.astype(np.float32)
    z = f(X, Y)
    # mmap_z[:] = f(mmap_xx, mmap_yy)

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
    # large_range_span = False
    try:
        z[np.abs(z) > z_val] = np.nan
        # mmap_z[np.abs(mmap_z) > z_val] = np.nan
        # Z = np.round(z, 2)  # Adjust precision as necessary
        # Replace inf with nan to avoid issues in contouring
        z[np.isinf(z)] = np.nan
        # CS = plt.contour(X, Y, np.ma.masked_where(z>z_val,
        #     z), levels=[0])
        # CS = plt.contour(X, Y, np.ma.masked_invalid(
        #     z), levels=[0], colors='blue', alpha=0)

        # z = adaptive_patch(f, z, 10, num_x, num_y, _x, _y)
        # z = adaptive_patch(f, z, 4)

        # z_masked = np.ma.masked_where(z < -z_val, z)

        # z[np.abs(z) < 1e-4] = 0.0

        # cont_gen = contour_generator(
        #     X, Y, z, name="serial")

        # 3. Populate with large data (example)
        # In reality, you'd fill this from your data source
        # mmap_z[:] = z
        cont_gen = contour_generator(
            X, Y, z, name="serial")
        # cont_gen = contour_generator(X, Y, z, quad_as_tri=True, corner_as_point=True, name="serial")
        del z
        # del z_masked
        del X
        del Y
        del _x
        del _y

        # lines(level) returns a list of branches (each is an (N, 2) array of coordinates)
        lines = cont_gen.lines(0)

        # 5. Clean up
        # del mmap_z  # Closes the file
        # del mmap_x  # Closes the file
        # del mmap_y  # Closes the file
        # del mmap_xx  # Closes the file
        # del mmap_yy  # Closes the file

        del cont_gen
        cont_gen = None

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
                    valid_mask = None
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
                valid_mask = None
                continue
            # if not large_range_span:
            #     max_y = np.max(segment[:, 1])
            #     min_y = np.min(segment[:, 1])
                # if np.abs(max_y - min_y) > 1e16:
                #     large_range_span = True

            all_points.append(base64.b64encode(
                segment.tobytes()).decode('utf-8'))
            # all_points.append(base64.b64encode(
            #     segment.astype(np.float32).tobytes()).decode('utf-8'))
            valid_mask = None
            del segment

        del lines

        # Check gc.garbage for uncollectable objects
        # print(f"Uncollectable objects: {gc.garbage}")

        # for obj in gc.garbage:
        #     print(f"Leaking object: {obj}, Type: {type(obj)}")
        # Optional: further inspect referrers
        # print(f"Referrers: {gc.get_referrers(obj)}")

        # Clear the garbage list for future collections if needed
        # gc.garbage.clear()

        # # # Disable debugging
        # gc.set_debug(0)

        return all_points, has_discontinuity, None

    except Exception as e:
        print(f"Error generating implicit plot points: {e}")
        return []
