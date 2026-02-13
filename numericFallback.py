from degree_radian import sin_mode, cos_mode, tan_mode, cot_mode, sec_mode, csc_mode, asin_mode, acos_mode, atan_mode, acot_mode, asec_mode, acsc_mode, trig_substitutions
from domain_finder import closer_boundary
import matplotlib.pyplot as plt
import sympy as sp
from sympy import symbols, solve, plot_implicit
import numpy as np
from scipy.optimize import fsolve, root, root_scalar

import matplotlib
matplotlib.use('Agg')


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


# Define variables
x, y = symbols('x y')


def generate_implicit_plot_points(expr, x_min=-10.0, x_max=10.0, y_min=-10.0, y_max=10.0,
                                  resolution=40000, adaptive=False, remove_temp_file=True):

    # Estimate y-range from symbolic solutions if possible
    try:
        Fy = sp.diff(expr, y)
        Fx = sp.diff(expr, x)
        critical_points = solve([Fy, Fx], (x, y))
        # print(f"Critical points: {critical_points}")
        if not critical_points.has(sp.I):
            if isinstance(critical_points, list):
                for cp in critical_points:
                    _e = expr.subs(x, cp[0])
                    res = solve(_e, y)
                    for r in res:
                        y_min = min(y_min, r)
                        y_max = max(y_max, r)
            else:
                _e = expr.subs(x, critical_points[x])
                res = solve(_e, y)
                for r in res:
                    y_min = min(y_min, float(r))
                    y_max = max(y_max, float(r))

    except:
        pass

    _x = np.linspace(x_min, x_max, 400)
    _y = np.linspace(y_min, y_max, 400)

    X, Y = np.meshgrid(_x, _y)
    # z = x**2 + y**2 - 1  # Example: circle equation x^2 + y^2 = 1
    f = sp.lambdify((x, y), expr, modules='numpy')
    z = f(X, Y)
    z[np.abs(z) > 40] = np.nan
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
                all_points.append(segment.tolist())
                # print("Extracted points for a contour segment:")
                # print(segment)

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
        if abs(slope) == abs(float('inf')) or abs(slope) > 15:
            infinite_discont = True
            if np.sign(branch[0][1]) == 1:
                branch.insert(0, [branch[0][0], "##"])
            else:
                branch.insert(0, [branch[0][0], "-##"])

        x0, y0 = branch[-2]
        x1, y1 = branch[-1]
        slope = (y1-y0)/(x1-x0) if x1 != x0 else float('inf')
        if abs(slope) == abs(float('inf')) or abs(slope) > 15:
            infinite_discont = True
            if np.sign(branch[-1][1]) == 1:
                branch.append([branch[-1][0], "##"])
            else:
                branch.append([branch[-1][0], "-##"])

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
