# %%fmt: off
import base64
from custom import custom
from sympy import factor, lambdify
from sympy.functions.elementary.trigonometric import TrigonometricFunction
import numpy as np
from sympy import symbols, solve
import sympy as sp
from my_misc import has_infinite_discontinuity_in_xrange, sanitize_contour_segments, generate_contoured_mesh, generate_contoured_mesh_fast
import re
import math
from contourpy import contour_generator

# Enable saving of uncollectable objects
# gc.set_debug(gc.DEBUG_SAVEALL)


# Define variables
# x, y = symbols('x y')
y = symbols('y')


def has_cusp(expr, _var, x_min=-10, x_max=10, y_min=-10, y_max=10):
    x = sp.Symbol(_var)
    if expr.has(TrigonometricFunction):
        return []
    try:
        pt = sp.solve([expr, sp.diff(expr, x), sp.diff(expr, y)],
                      [x, y], dict=True)
        if isinstance(pt, list):
            pt = [p for p in pt if x_min <= p[x] <=
                  x_max and y_min <= p[y] <= y_max]
        if len(pt) > 0:
            _v = expr.subs(x, pt[0][x]).subs(y, pt[0][y])
            if abs(_v) < 1e-12:
                return [pt[0][x], pt[0][y]]
        return []
    except Exception:
        return []


def find_cusp_points(expr, _var, x_min, x_max, y_min, y_max,
                     num_grid=250, tolerance=1e-6, max_results=50):
    """Return points (x, y) in the region that look like implicit curve cusps.

    A cusp candidate is a point on F(x,y)=0 where both partial derivatives
    vanish (Fx=0, Fy=0). The function attempts symbolic root finding first,
    then falls back to a grid-based numeric approximation.

    Args:
        expr: sympy expression for F(x,y)=0, or string parseable by sympy.
        _var: the variable name for the x-axis.
        x_min, x_max, y_min, y_max: search bounds for cusp points.
        num_grid: grid size for fallback numerical scan.
        tolerance: absolute tolerance for zero checks.
        max_results: max points to return.

    Returns:
        list of (x,y) tuples (floats) where cusps were found.
    """

    # res = cusp(expr)

    x = sp.Symbol(_var)

    if isinstance(expr, str):
        expr = sp.sympify(expr)

    if not isinstance(expr, sp.Expr):
        raise TypeError("expr must be a sympy expression or expression string")

    fx = sp.diff(expr, x)
    fy = sp.diff(expr, y)

    candidates = []

    # 1) Symbolic solve if available
    try:
        soln = sp.solve([expr, fx, fy], [x, y], dict=True)
        for item in soln:
            if x in item and y in item:
                px = float(item[x])
                py = float(item[y])
                if x_min <= px <= x_max and y_min <= py <= y_max:
                    candidates.append((px, py))
    except Exception:
        pass

    # Special handling for trigonometric functions to find all cusps
    if expr.has(TrigonometricFunction) and len(candidates) < max_results:
        try:
            f_at_y0 = expr.subs(y, 0)
            x_sols = sp.solveset(f_at_y0, x, sp.S.Reals)

            # Prepare lambdified functions for derivative validation
            fx_fun = lambdify((x, y), fx, modules=[custom, 'numpy'])
            fy_fun = lambdify((x, y), fy, modules=[custom, 'numpy'])

            if isinstance(x_sols, sp.ImageSet):
                lam = x_sols.lam
                base = lam.expr
                n = lam.variables[0]
                if sp.simplify(base / n) == sp.pi:
                    k_min = int(math.ceil(x_min / math.pi))
                    k_max = int(math.floor(x_max / math.pi))
                    for k in range(k_min, k_max + 1):
                        px = float(k * math.pi)
                        if x_min <= px <= x_max:
                            # Verify this is actually a cusp: both derivatives must be close to zero
                            try:
                                fx_val = abs(float(fx_fun(px, 0.0)))
                                fy_val = abs(float(fy_fun(px, 0.0)))
                                if fx_val < tolerance and fy_val < tolerance:
                                    candidates.append((px, 0.0))
                            except:
                                pass
            elif isinstance(x_sols, sp.Union) and len(x_sols.args) == 2 and all(isinstance(s, sp.ImageSet) for s in x_sols.args):
                # Special case for expressions like sin(x)**2 = 0, which has Union of even and odd multiples
                k_min = int(math.ceil(x_min / math.pi))
                k_max = int(math.floor(x_max / math.pi))
                for k in range(k_min, k_max + 1):
                    px = float(k * sp.pi)
                    if x_min <= px <= x_max:
                        # Verify this is actually a cusp: both derivatives must be close to zero
                        try:
                            fx_val = abs(float(fx_fun(px, 0.0)))
                            fy_val = abs(float(fy_fun(px, 0.0)))
                            if fx_val < tolerance and fy_val < tolerance:
                                candidates.append((px, 0.0))
                        except:
                            pass
            elif isinstance(x_sols, sp.FiniteSet):
                for sol in x_sols:
                    if sol.is_real and sol.is_finite:
                        px = float(sol)
                        if x_min <= px <= x_max:
                            # Verify this is actually a cusp: both derivatives must be close to zero
                            try:
                                fx_val = abs(float(fx_fun(px, 0.0)))
                                fy_val = abs(float(fy_fun(px, 0.0)))
                                if fx_val < tolerance and fy_val < tolerance:
                                    candidates.append((px, 0.0))
                            except:
                                pass
        except Exception:
            pass

    # 2) Numeric fallback (grid + approx)
    if len(candidates) < max_results and not expr.has(TrigonometricFunction):
        f_fun = lambdify((x, y), expr, modules=[custom, 'numpy'])
        fx_fun = lambdify((x, y), fx, modules=[custom, 'numpy'])
        fy_fun = lambdify((x, y), fy, modules=[custom, 'numpy'])

        xs = np.linspace(x_min, x_max, num_grid)
        ys = np.linspace(y_min, y_max, num_grid)
        X, Y = np.meshgrid(xs, ys, indexing='xy')

        with np.errstate(all='ignore'):
            F = f_fun(X, Y)
            Fx = fx_fun(X, Y)
            Fy = fy_fun(X, Y)

        # lambdify may return scalar when derivative is constant, so broadcast to grid shape
        if np.isscalar(F):
            F = np.full(X.shape, F, dtype=float)
        if np.isscalar(Fx):
            Fx = np.full(X.shape, Fx, dtype=float)
        if np.isscalar(Fy):
            Fy = np.full(X.shape, Fy, dtype=float)

        mask = np.isfinite(F) & np.isfinite(Fx) & np.isfinite(Fy)
        if np.any(mask):
            absF = np.abs(F[mask])
            absFx = np.abs(Fx[mask])
            absFy = np.abs(Fy[mask])

            # Skip if derivatives are never small enough for cusps
            if np.min(absFx) > 1e-3 or np.min(absFy) > 1e-3:
                pass  # No cusps possible
            else:
                # heuristics to choose likely cusp region
                f_th = max(tolerance, np.percentile(absF, 1) * 20)
                fx_th = max(tolerance, np.percentile(absFx, 1) * 20)
                fy_th = max(tolerance, np.percentile(absFy, 1) * 20)

                cand_mask = (
                    np.abs(F) <= f_th
                ) & (np.abs(Fx) <= fx_th) & (np.abs(Fy) <= fy_th)
                cand_mask &= mask

                indices = np.argwhere(cand_mask)
                for i, j in indices:
                    px = float(xs[j])
                    py = float(ys[i])
                    if x_min <= px <= x_max and y_min <= py <= y_max:
                        candidates.append((px, py))

    # dedupe close candidates
    uniq = []
    for p in candidates:
        if any(np.hypot(p[0] - q[0], p[1] - q[1]) < max(tolerance, 1e-5) for q in uniq):
            continue
        uniq.append(p)

    return uniq[:max_results]


def estimate_y_bounds2(equation, _var, x_min, x_max, num_x=400, y_min=None, y_max=None, y_samples=400, match_tol=None, f_tol=1e-15):
    # if (equation.has(TrigonometricFunction)) or ("_mode" in str(equation)):
    #     return (-300, 300)

    x = sp.Symbol(_var)

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


def grid_x_y_z_val(expr, _var, x_min, x_max, y_min, y_max):
    x = sp.Symbol(_var)

    z_val = 8
    has_discontinuity = has_infinite_discontinuity_in_xrange(
        expr, _var, x_min, x_max)

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


def generate_implicit_plot_points(expr, _var, x_min=-10.0, x_max=10.0, autoScale=False, has_discontinuity=False, y_min=-10.0, y_max=10.0):
    x = sp.Symbol(_var)

    # 2. Factor the expression to see the components
    factored_expr = sp.factor(expr)  # Result: x*(x - y)

    # Check if x is in the factors
    has_x = x in factored_expr.args

    if has_x:
        for arg in factored_expr.args:
            if arg.is_Number:
                continue
            elif arg == x:
                continue
            elif (arg.is_polynomial(y) and sp.degree(arg, gen=y) != 1) or (arg.is_polynomial(x) and sp.degree(arg, gen=x) != 1):
                continue
            else:
                expr = factored_expr/x

    has_1_divided_by_x = (1/x) in factored_expr.args

    if has_1_divided_by_x:
        for arg in factored_expr.args:
            if arg.is_Number:
                continue
            elif arg == 1/x:
                continue
            elif (arg.is_polynomial(y) and sp.degree(arg, gen=y) != 1) or (arg.is_polynomial(x) and sp.degree(arg, gen=x) != 1):
                continue
            else:
                expr = factored_expr*x

    (y_min, y_max) = estimate_y_bounds2(expr, _var, x_min, x_max)
    if y_min == 0:
        y_min = -10000
    if y_max == 0:
        y_max = 10000

    huge = False
    if abs(y_max - y_min) > 1e16:
        huge = True

    _max = np.max([np.abs(y_min), np.abs(y_max)])

    y_min = -_max
    y_max = _max

    cusp = find_cusp_points(expr, _var, x_min, x_max, y_min, y_max)
    # cusp = has_cusp(expr, _var, x_min, x_max, y_min, y_max)

    num_x, num_y, z_val, has_discontinuity = grid_x_y_z_val(
        expr, _var, x_min, x_max, y_min, y_max)
    if has_discontinuity:
        y_min = -500
        y_max = 500
    try:
        z_val = float(z_val)
    except Exception:
        pass

    f = lambdify((x, y), expr, modules=["numpy", custom])

    _x, _y = None, None
    deg_poly_y = 1
    if expr.is_polynomial(y):
        deg_poly_y = min([sp.degree(expr, gen=y), 50])

    deg_poly_x = 1
    if expr.is_polynomial(x):
        deg_poly_x = min([sp.degree(expr, gen=x), 50])

    # Pattern for common trig functions
    # pattern = r"\b(sin|cos|tan|asin|acos|atan|sin_mode|cos_mode|tan_mode|asin_mode|acos_mode|atan_mode)\b"

    # # Replace all matched trig functions with '1'
    # str_expr = re.sub(pattern, "1*", str(expr))
    # str_expr = sp.sympify(str_expr)

    # deg_poly_x = 1

    # if str_expr.is_polynomial(x):
    #     deg_poly_x = int(sp.degree(str_expr, x))

    # deg_poly_y = 1

    # if str_expr.is_polynomial(y):
    #     deg_poly_y = int(sp.degree(str_expr, y))

    num_x = 1000
    num_y = 1000

    if huge:
        _x = np.linspace(x_min, x_max, num_x)
        _f = 0.095
        _y = np.linspace(y_min, -1e12, int(num_y*_f), endpoint=False)
        _y = np.append(_y, np.linspace(-1e12, 1e12,
                                       int(num_y*(1-2*_f)), endpoint=False))
        _y = np.append(_y, np.linspace(1e12, y_max, int(num_y*_f)))
    elif deg_poly_y > 3:
        poly_x_factor = 1
        if deg_poly_x >= 2:
            poly_x_factor = 2.0
        _x = np.linspace(x_min, x_max, int(poly_x_factor*num_x*0.125))
        _y = np.linspace(y_min, y_max, num_y*8)
    else:
        initial_sz = 120
        try:
            _x, _y = generate_contoured_mesh_fast(
                f, [x_min, x_max], [y_min, y_max], initial_sz, 0.2)
        except Exception:
            _x, _y = generate_contoured_mesh(
                f, [x_min, x_max], [y_min, y_max], initial_sz, 0.2)

    density = len(_x)*len(_y)

    if density > 3000*3000:
        _x = np.linspace(x_min, x_max, 1600)
        _y = np.linspace(y_min, y_max, 900)

    if density < 400*400:
        _x = np.linspace(x_min, x_max, 1600)
        _y = np.linspace(y_min, y_max, 900)

    # if len(_x) > 3000:
    #     _x = np.linspace(x_min, x_max, 3000)
    # if len(_y) > 3000:
    #     _y = np.linspace(y_min, y_max, 3000)

    # if len(_x) < 400:
    #     _x = np.linspace(x_min, x_max, 400)
    # if len(_y) < 400:
    #     _y = np.linspace(y_min, y_max, 400)

    X, Y = np.meshgrid(_x, _y)
    del _x
    del _y

    z = f(X, Y)

    # z_val = 0.5

    if has_discontinuity:
        # z_val = float(2 * deg_poly_y)
        # z_val = float(2)
        z_val = 2.0
        # z_val = np.maximum(z_val, 2)

    try:
        # Ensure z is a numeric numpy array before masking values.
        z = np.asarray(z, dtype=np.float64)
        z[np.isinf(z)] = np.nan
        # Avoid invalid comparisons when z contains NaN values.
        with np.errstate(invalid='ignore'):
            mask = np.abs(z) > z_val
        z[mask] = np.nan

        cont_gen = contour_generator(
            X, Y, z, name="serial")
        # cont_gen = contour_generator(X, Y, z, quad_as_tri=True, corner_as_point=True, name="serial")
        del z
        # del z_masked
        del X
        del Y
        # del _x
        # del _y

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
                expr, _var, segment, x_min, x_max, has_discontinuity)
            if segment is None:
                valid_mask = None
                continue

            # cusp = []
            for _cusp in cusp:
                new_point = np.array([float(_cusp[0]), float(_cusp[1])])
                # Calculate distances from the new point to all points in the segment
                distances = np.sum((segment - new_point)**2, axis=1)
                # Find the index of the closest point in the segment to the new point
                closest_index = np.argmin(distances)
                # Insert the new point into the segment at the position of the closest point
                if closest_index == 0:
                    # _, idx = np.unique(segment, axis=0, return_index=True)
                    # segment = segment[np.sort(idx)]
                    segment = np.delete(segment, len(segment)-1, axis=0)
                    segment = np.delete(segment, 0, axis=0)
                    segment = np.insert(
                        segment, closest_index, new_point, axis=0)
                    segment = np.append(
                        segment, [new_point], axis=0)

                else:
                    # segment = np.insert(
                    #     segment, closest_index+1, new_point, axis=0)
                    # space1 = abs(segment[closest_index, 0] -
                    #              segment[closest_index-1, 0])
                    # space2 = abs(
                    #     segment[closest_index+1, 0] - segment[closest_index, 0])
                    # if abs(space2 > 4*space1):
                    #     continue

                    if abs(segment[closest_index, 1] - new_point[1]) > 1e-2:
                        # print(segment[closest_index, 0])
                        segment[closest_index] = new_point

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
