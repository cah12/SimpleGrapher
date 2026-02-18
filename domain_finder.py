import numpy as np
from scipy.optimize import brentq, newton
import warnings
import sympy as sp
from sympy import lambdify
from scipy.optimize import fsolve, root, root_scalar


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
def closer_boundary(fn, current_boundary, next_point, forward=True):
    y_step = np.abs(current_boundary[1] - next_point[1])
    x_step = np.abs(current_boundary[0] - next_point[0])

    # y_step = np.max([y_step, 0.01])
    # x_step = np.max([x_step, 0.01])

    stepFactor = 4

    y_range = (current_boundary[1]-stepFactor*y_step,
               current_boundary[1]+stepFactor*y_step)
    x_range = (current_boundary[0]-stepFactor*x_step,
               current_boundary[0]+1000*stepFactor*x_step)

    # if not forward:
    #     y_range = (current_boundary[1]+stepFactor*y_step,
    #                current_boundary[1]-stepFactor*y_step)
    #     x_range = (current_boundary[0]+stepFactor*x_step,
    #                current_boundary[0]-stepFactor*x_step)
    # elif forward:
    #     y_range = (current_boundary[0]-6*y_step, current_boundary[1]+6*y_step)
    #     x_range = (current_boundary[0]-6*x_step, current_boundary[1]+6*x_step)

    x_samples = np.linspace(x_range[0], x_range[1], 50)
    # y_samples = np.linspace(y_range[0], y_range[1], 1000)
    y_samples = np.linspace(-1e+4, 1e+4, 100)


   
    for x in x_samples:
        # print(x, is_x_in_domain_numerical(circle_func, x, y_range))
        if is_x_in_domain_numerical(fn, x, y_range):
            y = sp.symbols('y')
            _fn = fn(x, y)

            for i in range(len(y_samples) - 1):
                if _fn.subs(y, y_samples[i]) * _fn.subs(y, y_samples[i+1]) < 0:
                    return (x.item(), y_samples[i].item())

    return None


def generate_points_all_branches(f, x_min, x_max, num_x=400, y_min=None, y_max=None, y_samples=400, match_tol=None, f_tol=1e-15):
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
    x, y = sp.symbols('x y')
    x_vals = np.linspace(x_min, x_max, num_x)
    # f = lambdify((x, y), equation, 'numpy')

    # Estimate y-range from symbolic solutions if possible
    if y_min is None or y_max is None:
        try:
            y_sols = sp.solve(equation, y)
        except Exception:
            y_sols = []

        if y_sols:
            y_vals_est = []
            for ys in y_sols:
                try:
                    y_fun = lambdify(x, ys, 'numpy')
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

    return branches


# def closer_boundary(fn, current_boundary, next_point, forward=True):
#     y_step = np.abs(current_boundary[1] - next_point[1])
#     x_step = np.abs(current_boundary[0] - next_point[0])

#     # y_step = np.max([y_step, 0.01])
#     # x_step = np.max([x_step, 0.01])

#     stepFactor = 4

#     y_range = (current_boundary[0]-stepFactor*y_step,
#                current_boundary[1]+stepFactor*y_step)
#     x_range = (current_boundary[0]-stepFactor*x_step,
#                current_boundary[0]+stepFactor*x_step)

#     if not forward:
#         y_range = (current_boundary[0]+stepFactor*y_step,
#                    current_boundary[1]-stepFactor*y_step)
#         x_range = (current_boundary[0]+stepFactor*x_step,
#                    current_boundary[1]-stepFactor*x_step)
#     branches = generate_points_all_branches(fn, x_range[0], x_range[1])
#     return None


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
