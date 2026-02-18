from my_misc import estimate_y_bounds, sanitize_contour_segments
from degree_radian import sin_mode, cos_mode, tan_mode, cot_mode, sec_mode, csc_mode, asin_mode, acos_mode, atan_mode, acot_mode, asec_mode, acsc_mode, trig_substitutions
from domain_finder import closer_boundary
import matplotlib.pyplot as plt
import sympy as sp
from sympy import symbols, solve, plot_implicit
import numpy as np
from scipy.optimize import fsolve, root, root_scalar
from sympy.functions.elementary.trigonometric import TrigonometricFunction

from sympy import lambdify

from custom import custom

import matplotlib
matplotlib.use('Agg')



# Define variables
x, y = symbols('x y')


def estimate_y_bounds2(equation, x_min, x_max, num_x=400, y_min=None, y_max=None, y_samples=400, match_tol=None, f_tol=1e-15):
    if (equation.has(TrigonometricFunction)) or ("_mode" in str(equation)):
        return (-100, 100)
    
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

    _x = np.linspace(x_min, x_max, 1000)
    _y = np.linspace(y_min, y_max, 1000)

    X, Y = np.meshgrid(_x, _y)
    # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
    f = sp.lambdify((x, y), expr, modules=[custom, 'numpy'])
    z = f(X, Y)
    z_val = 0.5*y_max
    # Convert the expression to a string
    expr_str = str(expr)
    if (expr.has(TrigonometricFunction)) or ("_mode" in expr_str):
        z_val = np.maximum(z_val,100)
    else:
        z_val = np.maximum(z_val, 15)
    # if abs(z_val) >= 1e100:
    #     z_val = 100
    z[np.abs(z) > z_val] = np.nan
    try:
        # Z = np.round(z, 2)  # Adjust precision as necessary
        # Replace inf with nan to avoid issues in contouring
        z[np.isinf(z)] = np.nan
        CS = plt.contour(X, Y, np.ma.masked_invalid(
            z), levels=[0], colors='blue', alpha=0)

        all_points = []
        for level_segments in CS.allsegs:
            for segment in level_segments:
                # segment is a NumPy array of shape (n_points, 2), where each row is [x, y]
                all_points.append(segment)
                # all_points.append(segment.tolist())
                

        all_points = sanitize_contour_segments(expr, all_points, x_min, x_max)
        return all_points

    except Exception as e:
        print(f"Error generating implicit plot points: {e}")
        return []


def split_list(lst, delimiter):
    if delimiter not in lst:
        return [lst]
    else:
        return [lst[:lst.index(delimiter)], *split_list(lst[lst.index(delimiter)+1:], delimiter)]


def estimateStepSize(branch):
    previous = 1e300

    for i in range(1, len(branch)):
        x0 = branch[i-1][0]
        x1 = branch[i][0]
        if x0 == "##" or x1 == "##" or x0 == "#" or x1 == "#":
            continue
        s = x1 - x0
        if s < previous and s != 0:
            previous = abs(x1 - x0)
    return previous


def changeInSlope(branch, ind):

    x0 = branch[ind-1][0]
    x1 = branch[ind][0]
    y0 = branch[ind-1][1]
    y1 = branch[ind][1]

    slope1 = (y1-y0)/(x1-x0)
    ind = ind+1
    x0 = branch[ind-1][0]
    x1 = branch[ind][0]
    y0 = branch[ind-1][1]
    y1 = branch[ind][1]
    slope2 = (y1-y0)/(x1-x0)

    return slope2 - slope1


def points_with_vertical_tangent(f):
    # 2. Calculate partial derivatives
    df_dx = sp.diff(f, x)
    df_dy = sp.diff(f, y)

    # 3. Find where F_y = 0 and F = 0
    # For a circle, y=0 gives vertical tangents
    vertical_tangent_points = sp.solve([f, df_dy], (x, y))
    print(f"Points with vertical tangents: {vertical_tangent_points}")

    # 4. Numerical evaluation
    numerical_points = [(p[0].evalf(), p[1].evalf())
                        for p in vertical_tangent_points]
    print(f"Numerical points: {numerical_points}")

    return numerical_points


def verified(f, x_val, sign=1):
    # return True
    if f.has(TrigonometricFunction):
        return True
    try:
        v = f.subs({x: x_val, y: sign*1e+100})
        if abs(v) <= 1:
            return True
        return False
    except Exception:
        return False


def processBranches(branches, f):
    # return branches
    processed = []

    infinite_discont = False

    for branch in branches:
        if len(branch) < 40:
            continue

        # Find and remove straight lines
        # x0, y0 = branch[0]
        # x1, y1 = branch[-1]
        # ref_slope = (y1-y0)/(x1-x0) if x1 != x0 else float('inf')
        # a_straight_line = True
        # for i in range(1, len(branch)):
        #     x0, y0 = branch[i-1]
        #     x1, y1 = branch[i]
        #     slope = (y1-y0)/(x1-x0) if x1 != x0 else float('inf')
        #     # v = slope - ref_slope
        #     if abs(slope - ref_slope) > 0.01:
        #         a_straight_line = False
        #         break
        # if a_straight_line:
        #     continue

        # Mark infinity with "##"
        # vt = points_with_vertical_tangent(f)

        x0, y0 = branch[0]
        x1, y1 = branch[1]
        slope = (y1-y0)/(x1-x0) if x1 != x0 else float('inf')
        if abs(slope) == abs(float('inf')) or abs(slope) > 30:
            if np.sign(branch[0][1]) == 1 and verified(f, branch[0][0], 1):
                branch.insert(0, [branch[0][0], "##"])
                infinite_discont = True
            elif np.sign(branch[0][1]) == -1 and verified(f, branch[0][0], -1):
                branch.insert(0, [branch[0][0], "-##"])
                infinite_discont = True

        x0, y0 = branch[-2]
        x1, y1 = branch[-1]
        slope = (y1-y0)/(x1-x0) if x1 != x0 else float('inf')
        if abs(slope) == abs(float('inf')) or abs(slope) > 30:
            if np.sign(branch[-1][1]) == 1 and verified(f, branch[-1][0], 1):
                branch.append([branch[-1][0], "##"])
                infinite_discont = True
            elif np.sign(branch[-1][1]) == -1 and verified(f, branch[-1][0], -1):
                branch.append([branch[-1][0], "-##"])
                infinite_discont = True

        processed.append(branch)

    return processed, infinite_discont

##############################################################


def find_discontinuities_in_branch(branch, jump_factor=50, large_y_mult=1e3):
    """
    Find discontinuities in a single branch returned by `generate_points_all_branches`.

    Args:
        branch: list of [x, y] points (x increasing)
        jump_factor: multiplier for median dy to decide jump threshold
        large_y_mult: multiplier for median|y| to decide "infinite" threshold

    Returns:
        List of discontinuities where each entry is [index_in_branch_where_discontinuity_occur, type]
        - index is 1-based index of the earlier point involved in the discontinuity
        - type is one of: "infinite", "jump", "gap", "unknow"
    """
    if not branch or len(branch) < 2:
        return []

    ys = [pt[1] for pt in branch]
    try:
        ys_arr = np.asarray(ys, dtype=float)
    except Exception:
        return [[1, "unknow"]]

    dys = np.abs(np.diff(ys_arr))

    median_dy = float(np.median(dys)) if dys.size > 0 else 0.0
    median_abs_y = float(np.median(np.abs(ys_arr))) if ys_arr.size > 0 else 0.0

    jump_threshold = max(1e-12, median_dy * jump_factor)
    large_y_threshold = max(1e6, median_abs_y * large_y_mult)

    discontinuities = []

    for i in range(len(ys) - 1):
        y0 = ys_arr[i]
        y1 = ys_arr[i + 1]

        # gap: NaN present
        if np.isnan(y0) or np.isnan(y1):
            discontinuities.append([i + 1, "gap"])  # 1-based index
            continue

        # infinite-like: either value is infinite or very large compared to threshold
        if np.isinf(y0) or np.isinf(y1) or abs(y0) > large_y_threshold or abs(y1) > large_y_threshold:
            if abs(y1 - y0) > jump_threshold:
                discontinuities.append([i + 1, "infinite"])
            else:
                discontinuities.append([i + 1, "unknow"])
            continue

        # jump discontinuity
        if abs(y1 - y0) > jump_threshold:
            discontinuities.append([i + 1, "jump"])  # 1-based index
            continue

        # otherwise continuous at this segment

    return discontinuities
# Example: equation with y and x on both sides
# e.g., y**2 + x*y - x**2 = 0
# equation = y**2 + x*y - x**2
# equation = y - x**2
# equation = y**7+5*y**6-3*y**3+45 - x**2
# equation = y**7+5*y+45 - x**2
# equation = y**2 - x
# equation = y-1/sp.sin(x)


# Usage example:
# Symbolic approach
# points_sym = generate_points_symbolic(equation, x, y, -5, 5, 100)
# print("Symbolic solutions:")
# for i, branch in enumerate(points_sym):
#     print(f"Branch {i}: {len(branch)} points")

# Numerical approach using scipy.optimize.root
# points_num = generate_points_numerical(
#     equation, -100, 100, y_init_guess=1, num_points=1000, method='hybr')
# print(
#     f"\nNumerical solution (using scipy.optimize.root): {len(points_num)} points")
# print("First 5 points (list of lists format):")
# for point in points_num[:5]:
#     print(f"  {point}")

# Demo: y**2 - x = 0 -> two branches for x >= 0

# branches = generate_points_all_branches(
#     equation, -1, 10, num_x=300, y_samples=800)
# print(f"\nImplicit solver found {len(branches)} branches for y**2 - x = 0")
# for i, br in enumerate(branches):
#     print(f"Branch {i}: {len(br)} points, first 5: {br[:5]}")

# Demo on the earlier high-degree equation (smaller grid to keep runtime reasonable)
# branches2 = generate_points_all_branches(
#     equation, -10, 10, num_x=200, y_samples=600)
# print(
#     f"\nHigh-degree equation -> {len(branches2)} branches (sample up to first 3 shown)")
# for i, br in enumerate(branches2[:3]):
#     print(f"Branch {i}: {len(br)} points, sample first 5: {br[:5]}")
# print(
#     f"\nHigh-degree equation -> Number of branches: {len(branches2)}")
# for i, br in enumerate(branches2):
#     print(f"*********Branch {i}: {len(br)} points************")
#     for pt in br:
#         print(pt)
#     print("")

# points = merge_branches_with_discontinuities(branches2)

# for p in points:
#     print(p)
