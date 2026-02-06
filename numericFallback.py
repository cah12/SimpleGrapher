import sympy as sp
from sympy import symbols, solve, lambdify
import numpy as np
from scipy.optimize import fsolve, root, root_scalar

from domain_finder import closer_boundary

# Define variables
x, y = symbols('x y')


# Method 1: Using SymPy's solve (works for equations solvable for y)
def generate_points_symbolic(equation, x_var, y_var, x_min, x_max, num_points=100):
    """Solve equation for y symbolically, then evaluate numerically"""
    # Solve for y in terms of x
    y_solutions = solve(equation, y_var)

    x_vals = np.linspace(x_min, x_max, num_points)
    points = []

    for y_sol in y_solutions:
        y_func = lambdify(x_var, y_sol, 'numpy')
        try:
            y_vals = y_func(x_vals)
            # Filter out complex/NaN values
            valid = ~np.isnan(y_vals) & np.isreal(y_vals)
            points.append(list(zip(x_vals[valid], np.real(y_vals[valid]))))
        except:
            pass

    return points

# Method 2: Using scipy.optimize.root for implicit equations


def generate_points_numerical(equation, x_min, x_max, y_init_guess=0, num_points=100, method='hybr'):
    """
    Generate x-y points for implicit equations using scipy.optimize.root.

    Args:
        equation: SymPy equation (implicitly = 0)
        x_min: Lower bound for x
        x_max: Upper bound for x
        y_init_guess: Initial guess for y (default: 0)
        num_points: Number of points to generate (default: 100)
        method: Root finding method - 'hybr', 'lm', 'broyden1', etc. (default: 'hybr')

    Returns:
        List of lists [[x1, y1], [x2, y2], ...] representing points on the curve
    """
    x_vals = np.linspace(x_min, x_max, num_points)
    points = []

    # Create a function from the equation
    f = lambdify((x, y), equation, 'numpy')

    current_y = y_init_guess
    for x_val in x_vals:
        try:
            # Define function to find root: f(x_val, y) = 0
            def equation_at_x(y_val):
                return f(x_val, y_val)

            # Use scipy.optimize.root
            result = root(equation_at_x, current_y, method=method)

            if result.success:
                y_sol = result.x[0] if hasattr(
                    result.x, '__len__') else result.x
                points.append([x_val.item(), y_sol.item()])
                current_y = y_sol  # Use previous solution as next guess
        except:
            pass

    return points


def generate_points_all_branches(equation, x_min, x_max, num_x=400, y_min=None, y_max=None, y_samples=400, match_tol=None, f_tol=1e-15):
    """
    Robustly compute all real y(x) branches for the implicit equation `equation(x,y)=0` over x in [x_min, x_max].

    Strategy:
    - Try to use symbolic solutions to estimate y-range (if available).
    - For each x sample, evaluate f(x, y) on a dense y grid, find sign changes and near-zero samples to bracket roots.
    - Use `scipy.optimize.root_scalar` (brentq) to find roots for each bracket.
    - Group roots across successive x values into continuous branches.

    Returns:
        List of branches, where each branch is a list of [x, y] points (both floats).
    """
    x_vals = np.linspace(x_min, x_max, num_x)
    f = lambdify((x, y), equation, 'numpy')

    # Estimate y-range from symbolic solutions if possible
    # if y_min is None or y_max is None:
    #     try:
    #         y_sols = solve(equation, y)
    #     except Exception:
    #         y_sols = []

    #     if y_sols:
    #         y_vals_est = []
    #         for ys in y_sols:
    #             try:
    #                 y_fun = lambdify(x, ys, 'numpy')
    #                 y_eval = y_fun(x_vals)
    #                 y_eval = np.asarray(y_eval)
    #                 valid = ~np.isnan(y_eval) & np.isreal(y_eval)
    #                 if np.any(valid):
    #                     y_vals_est.append(np.real(y_eval[valid]))
    #             except Exception:
    #                 pass
    #         if y_vals_est:
    #             all_est = np.hstack(y_vals_est)
    #             est_min, est_max = float(
    #                 np.min(all_est)), float(np.max(all_est))
    #             padding = max(1.0, 0.1 * (est_max - est_min))
    #             if y_min is None:
    #                 y_min = est_min - padding
    #             if y_max is None:
    #                 y_max = est_max + padding

    # Fallback heuristic if still not set
    if y_min is None or y_max is None:
        y_guess = max(1.0, abs(x_min), abs(x_max)) * 10.0
        y_min = -y_guess if y_min is None else y_min
        y_max = y_guess if y_max is None else y_max

    # matching tolerance for grouping branches
    if match_tol is None:
        match_tol = (y_max - y_min) * 0.1 if (y_max - y_min) > 0 else 0.1

    branches = []  # list of lists of [x,y]
    prev_roots = []  # stores last y of each branch

    y_grid = np.linspace(y_min, y_max, y_samples)

    for xi in x_vals:
        try:
            fvals = f(xi, y_grid)
            fvals = np.asarray(fvals, dtype=float)
        except Exception:
            # If evaluation fails, skip this x
            continue

        roots_this_x = []
        # points where function is (close to) zero on grid
        close_idx = np.where(np.abs(fvals) <= f_tol)[0]
        for idx in close_idx:
            roots_this_x.append(float(y_grid[idx]))

        # bracketing sign changes
        sign_changes = np.where(
            np.sign(fvals[:-1]) * np.sign(fvals[1:]) < 0)[0]
        for idx in sign_changes:
            a, b = float(y_grid[idx]), float(y_grid[idx + 1])
            try:
                sol = root_scalar(lambda yy: float(f(xi, yy)),
                                  bracket=[a, b], method='brentq')
                if sol.converged:
                    roots_this_x.append(float(sol.root))
            except Exception:
                pass

        # unique roots (within tolerance)
        if not roots_this_x:
            # no real roots at this x
            prev_roots = []
            continue

        roots = sorted(set(np.round(np.array(roots_this_x, dtype=float), 12)))

        # match roots to existing branches (by proximity)
        assigned = [False] * len(roots)
        new_prev_roots = []

        # First, try to match to existing branches
        for bi, brow in enumerate(prev_roots):
            # find closest root to brow
            diffs = [abs(r - brow) for r in roots]
            if diffs:
                best_idx = int(np.argmin(diffs))
                if diffs[best_idx] <= match_tol and not assigned[best_idx]:
                    # append to branch bi
                    branches[bi].append([float(xi), float(roots[best_idx])])
                    assigned[best_idx] = True
                    new_prev_roots.append(roots[best_idx])
                else:
                    # branch disappears at this x (no append)
                    new_prev_roots.append(brow)
            else:
                new_prev_roots.append(brow)

        # Create new branches for unassigned roots
        for ri, r in enumerate(roots):
            if not assigned[ri]:
                branches.append([[float(xi), float(r)]])
                new_prev_roots.append(r)

        # Update prev_roots for next iteration
        prev_roots = new_prev_roots

    # Filter branches: keep only branches with at least 2 points
    branches = [br for br in branches if len(br) >= 2]

    # def fnc(x, y):
    #     return y**6+8*y-x
    # return equation  # return symbolic equation

    fnc = lambdify((x, y), equation, 'numpy')

    # current_boundary = (-6.784, -0.9276)
    # next_point = (-6.683, -0.9034)

    for branch in branches:
        if branch[0][0] != x_min:
            current_boundary = (branch[0][0], branch[0][1])
            next_point = (branch[1][0], branch[1][1])
            b = closer_boundary(
                fnc, current_boundary, next_point, True)
            if b is not None:
                branch.insert(0, [b[0], b[1]])

        if branch[len(branch) - 1][0] != x_max:
            current_boundary = (
                branch[len(branch) - 1][0], branch[len(branch) - 1][1])
            next_point = (branch[len(branch) - 2][0],
                          branch[len(branch) - 2][1])
            b = closer_boundary(
                fnc, current_boundary, next_point, False)
            if b is not None:
                branch.append([b[0], b[1]])

    return branches


def merge_branches_with_discontinuities(branches, gap_marker="~~", jump_marker="~", pos_inf="+oo", neg_inf="-oo", jump_factor=50, large_y_mult=1e3):
    """
    Merge branches into a single ordered list of points suitable for plotting, inserting discontinuity markers.

    Input:
        branches: list of branches, each branch is a list of [x (float), y (float)]. Branch x-values must be increasing.
    Output:
        List of [x, y] where y is either a float or a marker string ("+oo", "-oo", "~", "~~").

    Heuristics used:
    - Compute median step in y across all branches to set a jump threshold: if a step between successive points
      exceeds `median_dy * jump_factor`, a jump is declared.
    - If one of the values around a jump is very large (abs(y) > large_y_threshold), treat as infinite (+oo/-oo).
    - Otherwise treat as a finite jump (~). Missing or NaN values generate a gap marker (~~).
    - Branches are concatenated in ascending starting-x order and separated by a gap marker to avoid connecting them.
    """
    if not branches:
        return []

    # collect statistics
    all_dys = []
    all_abs_y = []
    for br in branches:
        ys = [pt[1] for pt in br]
        all_abs_y.extend([abs(y) for y in ys])
        all_dys.extend([abs(ys[i + 1] - ys[i]) for i in range(len(ys) - 1)])

    median_dy = float(np.median(all_dys)) if all_dys else 0.0
    median_abs_y = float(np.median(all_abs_y)) if all_abs_y else 0.0

    jump_threshold = max(1e-12, median_dy * jump_factor)
    large_y_threshold = max(1e6, median_abs_y * large_y_mult)

    def process_branch(br):
        out = []
        for i, (x, y) in enumerate(br):
            out.append([float(x), float(y)])
            if i + 1 < len(br):
                x2, y2 = br[i + 1]
                # handle NaNs
                if np.isnan(y) or np.isnan(y2):
                    midx = (x + x2) / 2.0
                    out.append([float(midx), gap_marker])
                    continue

                dy = abs(y2 - y)
                if dy > jump_threshold:
                    midx = (x + x2) / 2.0
                    # check for infinite-like behavior
                    if abs(y) > large_y_threshold or abs(y2) > large_y_threshold:
                        large_val = y if abs(y) >= abs(y2) else y2
                        marker = pos_inf if large_val > 0 else neg_inf
                        out.append([float(midx), marker])
                    else:
                        out.append([float(midx), jump_marker])
        return out

    processed = [process_branch(br) for br in branches if br]

    # sort branches by their starting x
    processed = sorted(processed, key=lambda br: br[0][0])

    merged = []
    for br in processed:
        if not br:
            continue
        if merged:
            # insert a gap marker between branches to avoid connecting them
            last_x = merged[-1][0]
            start_x = br[0][0]
            midx = (last_x + start_x) / 2.0
            merged.append([float(midx), gap_marker])
        merged.extend(br)

    # compress consecutive identical markers
    out = []
    prev_marker = None
    for x, y in merged:
        if isinstance(y, str):
            if y == prev_marker:
                continue
            prev_marker = y
            out.append([x, y])
        else:
            prev_marker = None
            out.append([x, float(y)])

    return out


##############################################################
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
