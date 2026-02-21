import gc
from my_misc import estimate_y_bounds, sanitize_contour_segments
from degree_radian import sin_mode, cos_mode, tan_mode, cot_mode, sec_mode, csc_mode, asin_mode, acos_mode, atan_mode, acot_mode, asec_mode, acsc_mode, trig_substitutions
from domain_finder import closer_boundary

import sympy as sp
from sympy import symbols, solve, plot_implicit
import numpy as np
from scipy.optimize import fsolve, root, root_scalar
from sympy.functions.elementary.trigonometric import TrigonometricFunction

from sympy import lambdify

from custom import custom

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



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


def generate_implicit_plot_points(expr, x_min=-10.0, x_max=10.0, y_min=-10.0, y_max=10.0,
                                  resolution=40000, adaptive=False, remove_temp_file=True):  
    
    


    (_y_min, _y_max) = estimate_y_bounds2(expr, x_min, x_max)
    # (_y_min, _y_max) = (x_min, x_max)


    y_min = min(y_min, _y_min)
    y_max = max(y_max, _y_max)

    _x = np.linspace(x_min, x_max, 350, dtype=np.float32)
    _y = np.linspace(y_min, y_max, 700, dtype=np.float32)

    X, Y = np.meshgrid(_x, _y)
    # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
    f = sp.lambdify((x, y), expr, modules=[custom, 'numpy'])
    z = f(X, Y)
    z_val = 0.03*y_max
    # Convert the expression to a string
    expr_str = str(expr)
    if expr.has(TrigonometricFunction):
        z_val = np.maximum(z_val,10)
    if ("_mode" in expr_str):
        z_val = np.maximum(z_val,13)
    else:
        z_val = np.maximum(z_val, 15)
    # if abs(z_val) >= 1e100:
    # z_val = 10
    z[np.abs(z) > z_val] = np.nan
    try:
        # Z = np.round(z, 2)  # Adjust precision as necessary
        # Replace inf with nan to avoid issues in contouring
        z[np.isinf(z)] = np.nan
        CS = plt.contour(X, Y, np.ma.masked_invalid(
            z), levels=[0])  # Use _legacy_new_subsegments=False for use_legacy_contour=True)
        # CS = plt.contour(X, Y, np.ma.masked_invalid(
        #     z), levels=[0], colors='blue', alpha=0)

        all_points = []
        all_segments = CS.allsegs  # Get the list of segments for the first contour level
        for level_segments in all_segments:
            for segment in level_segments:
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


                all_points.append(segment.tolist())
                segment = None  # Clear reference to segment to free memory
        all_segments = None  # Clear reference to all_segments to free memory    

        all_points = sanitize_contour_segments(expr, all_points, x_min, x_max)
        # del CS
        plt.clf()
        plt.cla()
        # gc.collect()
        # plt.close()
        # Mandatory cleanup
        
        plt.close('all')
        del CS
        gc.collect() # Force garbage collection
        return all_points

    except Exception as e:
        print(f"Error generating implicit plot points: {e}")
        return []

