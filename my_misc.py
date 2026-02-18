
import numpy as np
from typing import List, Tuple
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from domain_finder import is_x_in_domain_numerical, closer_boundary
from custom import custom
from sympy.functions.elementary.trigonometric import TrigonometricFunction


def sanitize_contour_segments(expr, allsegs: List[np.ndarray],
                               x_min: float = -1e6, x_max: float = 1e6,threshold_distance: float = 1e-6,
                              max_segment_length: float = 1e6) -> List[np.ndarray]:
    """
    Sanitize contour segments by removing spurious and rogue lines.

    Handles complex implicit functions f(x, y)=0, trigonometric discontinuities,
    and marks infinity with "##" at y-values.

    Parameters
    ----------
    allsegs : List[np.ndarray]
        Contour segments from contour level 0, where each segment is an (N, 2) array
        of (x, y) coordinates.
    threshold_distance : float, optional
        Minimum distance for a segment to be considered valid (default: 1e-6).
    max_segment_length : float, optional
        Maximum allowed segment length before flagging as potentially spurious (default: 1e6).

    Returns
    -------
    List[np.ndarray] or List[Tuple/List]
        Sanitized segments with infinity marked as "##" in y-value where appropriate.
        Each segment is either a numpy array of valid coordinates or contains "##" markers.
    """

    sanitized = []

    for segment in allsegs:
        if segment is None or len(segment) < 40:
            continue

        # Remove NaN values
        valid_mask = ~np.isnan(segment).any(axis=1)
        if not np.any(valid_mask):
            continue

        segment = segment[valid_mask]

        if len(segment) < 2:
            continue

        # Check for degenerate segments (all points identical or too close)
        distances = np.sqrt(np.sum(np.diff(segment, axis=0)**2, axis=1))
        if np.all(distances < threshold_distance):
            continue

        # Handle trigonometric discontinuities and large jumps
        cleaned_segment = _handle_discontinuities(expr, segment, x_min, x_max, max_segment_length)

        if cleaned_segment is not None:
            sanitized.append(cleaned_segment.tolist())

    return sanitized


def _handle_discontinuities(expr, segment: np.ndarray, x_min: float = -1e6, x_max: float = 1e6,
                            max_segment_length: float = 1e6) -> np.ndarray:
    """
    Handle discontinuities from trigonometric functions and infinity points.

    Splits segments at discontinuities where points are connected but should not be.
    Marks infinity with "##" in y-value.

    Parameters
    ----------
    segment : np.ndarray
        Segment with shape (N, 2) containing (x, y) pairs.
    max_segment_length : float
        Maximum allowed distance between consecutive points.

    Returns
    -------
    np.ndarray or List
        Cleaned segment, possibly with "##" markers for infinity.
    """
    

    if len(segment) < 2:
        return None
    
    if not has_infinite_discontinuity_in_xrange(expr, x_min, x_max):
        return segment

    # Calculate distances between consecutive points
    distances = np.sqrt(np.sum(np.diff(segment, axis=0)**2, axis=1))

    # Find discontinuities (large jumps)
    discontinuity_indices = np.where(distances > max_segment_length)[0]

    # If no major discontinuities, check for infinity markers
    if len(discontinuity_indices) == 0:
        return _mark_infinity_points(expr, segment)

    # Split segment at discontinuities
    split_segments = []
    start_idx = 0

    for disc_idx in discontinuity_indices:
        end_idx = disc_idx + 1

        # Extract sub-segment
        sub_segment = segment[start_idx:end_idx]

        if len(sub_segment) >= 2:
            # Check if the jump indicates an asymptote/infinity
            if _is_asymptote_jump(segment[disc_idx], segment[disc_idx + 1]):
                # Mark the endpoints with infinity indicator
                marked_sub = _mark_infinity_at_endpoints(expr, sub_segment)
                split_segments.append(marked_sub)
            else:
                split_segments.append(_mark_infinity_points(expr, sub_segment))

        start_idx = end_idx

    # Add remaining segment
    if start_idx < len(segment):
        sub_segment = segment[start_idx:]
        if len(sub_segment) >= 2:
            split_segments.append(_mark_infinity_points(expr, sub_segment))

    # Merge back if only one segment
    if len(split_segments) == 1:
        return split_segments[0]
    elif len(split_segments) > 1:
        return split_segments

    return _mark_infinity_points(segment)


def _is_asymptote_jump(point1: np.ndarray, point2: np.ndarray,
                       angle_threshold: float = 0.95) -> bool:
    """
    Detect if a large jump indicates a trigonometric or other asymptote.

    Parameters
    ----------
    point1, point2 : np.ndarray
        Consecutive points with potential asymptote between them.
    angle_threshold : float
        Cosine threshold for detecting near-vertical asymptotes (default: 0.95).

    Returns
    -------
    bool
        True if jump appears to be an asymptote.
    """

    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]

    # Check for near-infinite slope (vertical asymptote)
    if abs(dx) < 1e-8 and abs(dy) > 1:
        return True

    # Check for near-horizontal asymptote approach
    if abs(dy) > abs(dx) * 10:
        return True

    # Check if both x and y differences are huge (suggests function discontinuity)
    distance = np.sqrt(dx**2 + dy**2)
    if distance > 1e5:
        return True

    return False





def has_infinite_discontinuity_in_xrange(implicit_expr, x_min, x_max) -> bool:
    """
    Symbolically detect vertical/infinite discontinuities for f(x,y)=0
    between the x-values of `current_boundary` and `next_point`.

    Approach:
    - Convert `implicit_expr` to a SymPy expression `F(x,y)` (if possible).
    - Consider the numerator `N(x,y)` of `F` (so we inspect `N==0`).
    - If `N` is a polynomial in `y`, compute its leading coefficient in `y`.
      Solutions of `leading_coeff(x)==0` are x-values where the highest-power
      term in `y` vanishes; these are candidate x-locations where arbitrarily
      large |y| can satisfy `N==0` (vertical asymptote behavior).
    - As a fallback, compute `lim_{y->oo} F(x,y)` and solve for x where that
      limit equals 0.

    Returns True if any symbolic candidate lies within the x-range between
    the two provided boundary points; False otherwise.
    """
    # extract x-range from provided boundary points (expect sequences/tuples)
    def _get_x(val):
        try:
            return float(val[0])
        except Exception:
            try:
                return float(val)
            except Exception:
                raise ValueError("Cannot extract x coordinate from boundary values")

    x0 = x_min
    x1 = x_max
    xmin = min(x0, x1)
    xmax = max(x0, x1)

    # Prepare symbolic expression
    expr_sym = None
    if isinstance(implicit_expr, sp.Expr):
        expr_sym = implicit_expr
    elif isinstance(implicit_expr, str):
        s = implicit_expr.strip()
        if "=" in s:
            L, R = s.split("=", 1)
            s = f"({L})-({R})"
        try:
            expr_sym = parse_expr(s)
        except Exception:
            expr_sym = None

    x, y = sp.symbols("x y")
    if expr_sym is None:
        return False

    # Work with numerator (if rational) so we check polynomial condition N(x,y)==0
    try:
        N, D = sp.together(expr_sym).as_numer_denom()
    except Exception:
        N = expr_sym

    # Try polynomial-in-y route
    try:
        polyN = sp.Poly(sp.expand(N), y)
        deg = polyN.degree()
        if deg is not None and deg >= 1:
            leading = polyN.LC()
            # Solve leading(x) == 0 on reals
            solset = sp.solveset(sp.Eq(sp.simplify(sp.factor(leading)), 0), x, domain=sp.S.Reals)
            # If solution set intersects the interval, report True
            if isinstance(solset, sp.Set):
                interval = sp.Interval(float(xmin), float(xmax))
                inter = solset.intersect(interval)
                if not inter.is_EmptySet:
                    return True
            else:
                # fallback: iterate finite solutions
                for s in solset:
                    try:
                        sv = float(s.evalf())
                        if xmin - 1e-12 <= sv <= xmax + 1e-12:
                            return True
                    except Exception:
                        continue
    except Exception:
        # not polynomial in y or other failure; fall through to limit approach
        pass

    # Fallback: compute limit as y->oo of expr_sym; if limit simplifies to 0
    # for some x in the interval, that's indicative of vertical asymptote.
    try:
        lim_expr = sp.limit(expr_sym, y, sp.oo)
        # Solve lim_expr == 0 over reals
        solset2 = sp.solveset(sp.Eq(sp.simplify(sp.factor(lim_expr)), 0), x, domain=sp.S.Reals)
        if isinstance(solset2, sp.Set):
            interval = sp.Interval(float(xmin), float(xmax))
            inter2 = solset2.intersect(interval)
            if not inter2.is_EmptySet:
                return True
        else:
            for s in solset2:
                try:
                    sv = float(s.evalf())
                    if xmin - 1e-12 <= sv <= xmax + 1e-12:
                        return True
                except Exception:
                    continue
    except Exception:
        pass

    return False

# def _mark_infinity_points(segment: np.ndarray) -> np.ndarray:
#     """
#     Identify and mark points near infinity in a segment.

#     Parameters
#     ----------
#     segment : np.ndarray
#         Segment with shape (N, 2) containing (x, y) pairs.

#     Returns
#     -------
#     np.ndarray or List
#         Segment with "##" marking y-values of points at infinity.
#     """

#     # Use large numeric sentinels for infinity markers so arrays remain numeric
#     POS_INF = 1e300
#     NEG_INF = -1e300

#     result = []
#     for x, y in segment:
#         # Positive infinity or very large positive values
#         if np.isposinf(y) or (not np.isinf(y) and y > 1e15):
#             result.append([x, POS_INF])
#         # Negative infinity or very large negative values
#         elif np.isneginf(y) or (not np.isinf(y) and y < -1e15):
#             result.append([x, NEG_INF])
#         else:
#             result.append([x, y])

#     return np.array(result, dtype=float)

def _mark_infinity_points(expr, segment: np.ndarray) -> np.ndarray:
    """
    Identify and mark points near infinity in a segment.

    Parameters
    ----------
    segment : np.ndarray
        Segment with shape (N, 2) containing (x, y) pairs.

    Returns
    -------
    np.ndarray or List
        Segment with "##" marking y-values of points at infinity.
    """

    # Use large numeric sentinels for infinity markers so arrays remain numeric
    POS_INF = 1e300
    NEG_INF = -1e300
    THRESHOLD_SLOPE = 120

    if len(segment) < 2:
        return segment
    
    # _max_y = segment.max(axis=0)[1]
    # _min_y = segment.min(axis=0)[1]
    
    x_0, y_0 = segment[0]
    x_1, y_1 = segment[1]

    if (expr.has(TrigonometricFunction)) or ("_mode" in str(expr)):
        x_0 = np.deg2rad(x_0)
        x_1 = np.deg2rad(x_1)
    
    slope = (y_1 - y_0)/(x_1 - x_0)
    if abs(slope) > THRESHOLD_SLOPE:
        if np.sign(y_0) == -1:
            segment[0,1] = NEG_INF             
        else:
            segment[0,1] = POS_INF 
        # if np.sign(y_0) == -1 and y_1 > y_0:
        #     segment[0,1] = NEG_INF             
        # elif np.sign(y_0) == 1 and y_1 < y_0:
        #     segment[0,1] = POS_INF 
        
    x_0, y_0 = segment[len(segment)-1]
    x_1, y_1 = segment[len(segment)-2]
    if (expr.has(TrigonometricFunction)) or ("_mode" in str(expr)):
        x_0 = np.deg2rad(x_0)
        x_1 = np.deg2rad(x_1)
    slope = (y_1 - y_0)/(x_1 - x_0)
    if abs(slope) > THRESHOLD_SLOPE:
        if np.sign(y_0) == -1:
            segment[len(segment)-1,1] = NEG_INF 
        else:
            segment[len(segment)-1,1] = POS_INF 
        # if np.sign(y_0) == -1 and y_1 > y_0:
        #     segment[len(segment)-1,1] = NEG_INF 
        # elif np.sign(y_0) == 1 and y_1 < y_0:
        #     segment[len(segment)-1,1] = POS_INF 
    return segment        


# def _mark_infinity_at_endpoints(segment: np.ndarray) -> np.ndarray:
#     """
#     Mark the endpoints of a segment that appears to approach asymptotes.

#     Parameters
#     ----------
#     segment : np.ndarray
#         Segment with shape (N, 2).

#     Returns
#     -------
#     np.ndarray or List
#         Segment with endpoints marked for infinity where appropriate.
#     """

#     POS_INF = 1e300
#     NEG_INF = -1e300

#     result = []
#     for i, (x, y) in enumerate(segment):
#         # Mark last point if it's at the discontinuity
#         if i == len(segment) - 1 and (np.isposinf(y) or (not np.isinf(y) and y > 1e12)):
#             result.append([x, POS_INF])
#         elif i == len(segment) - 1 and (np.isneginf(y) or (not np.isinf(y) and y < -1e12)):
#             result.append([x, NEG_INF])
#         else:
#             result.append([x, y])

#     return np.array(result, dtype=float)


def _mark_infinity_at_endpoints(expr, segment: np.ndarray) -> np.ndarray:
    """
    Identify and mark points near infinity in a segment.

    Parameters
    ----------
    segment : np.ndarray
        Segment with shape (N, 2) containing (x, y) pairs.

    Returns
    -------
    np.ndarray or List
        Segment with "##" marking y-values of points at infinity.
    """

    # Use large numeric sentinels for infinity markers so arrays remain numeric
    POS_INF = 1e300
    NEG_INF = -1e300
    THRESHOLD_SLOPE = 300

    if len(segment) < 2:
        return segment
    
    x_0, y_0 = segment[0]
    x_1, y_1 = segment[1]
    new_row = None
    slope = (y_1 - y_0)/(x_1 - x_0)
    if abs(slope) > THRESHOLD_SLOPE:
        if np.sign(y_0) == -1:
            segment[0,1] = NEG_INF             
        else:
            segment[0,1] = POS_INF 
        
    x_0, y_0 = segment[len(segment)-1]
    x_1, y_1 = segment[len(segment)-2]
    slope = (y_1 - y_0)/(x_1 - x_0)
    if abs(slope) > THRESHOLD_SLOPE:
        if np.sign(y_0) == -1:
            segment[len(segment)-1,1] = NEG_INF 
        else:
            segment[len(segment)-1,1] = POS_INF 
    return segment  


# Example usage and testing
if __name__ == "__main__":
    # Example 1: Basic contour segments
    seg1 = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    seg2 = np.array([[0, 5], [1, 6], [2, 7]])

    allsegs = [seg1, seg2]
    sanitized = sanitize_contour_segments(allsegs)

    print("Sanitized segments:")
    for i, seg in enumerate(sanitized):
        print(f"Segment {i}:")
        print(seg)
        print()

    # Example 2: Segment with discontinuity (asymptote)
    seg_with_discontinuity = np.array(
        [[0, 0], [1, 1], [1.001, 1e10], [2, -1e10], [3, -1], [4, -2]])
    allsegs2 = [seg_with_discontinuity]
    sanitized2 = sanitize_contour_segments(allsegs2)

    print("Sanitized segments with discontinuity:")
    for i, seg in enumerate(sanitized2):
        print(f"Segment {i}:")
        print(seg)
        print()

    # Example 3: Segment with NaN values
    seg_with_nan = np.array([[0, 0], [1, np.nan], [2, 2], [3, 3]])
    allsegs3 = [seg_with_nan]
    sanitized3 = sanitize_contour_segments(allsegs3)

    print("Sanitized segments with NaN:")
    for i, seg in enumerate(sanitized3):
        print(f"Segment {i}:")
        print(seg)





def estimate_y_bounds(implicit_func, lower_x: float, upper_x: float,
                      num_samples: int = 50, y_search_range: Tuple[float, float] = (-1e100, 1e100)) -> Tuple[float, float]:
    """
    Estimate the Y bounds for an implicit function over a given X range.

    For an implicit function f(x, y) = 0, this function samples x values across
    [lower_x, upper_x] and finds the corresponding y values that satisfy the equation.

    Parameters
    ----------
    implicit_func : callable
        A function f(x, y) that should equal 0 at solution points. Signature: f(x, y) -> float
    lower_x : float
        The lower bound of the X range to search.
    upper_x : float
        The upper bound of the X range to search.
    num_samples : int, optional
        Number of x values to sample in the range (default: 50).
    y_search_range : Tuple[float, float], optional
        The range [min_y, max_y] to search for y solutions (default: (-1e6, 1e6)).

    Returns
    -------
    Tuple[float, float]
        A tuple (lower_y, upper_y) representing the estimated bounds of y values.

    Examples
    --------
    >>> # Circle: x^2 + y^2 - 1 = 0
    >>> circle = lambda x, y: x**2 + y**2 - 1
    >>> lower_y, upper_y = estimate_y_bounds(circle, -1, 1)
    >>> print(f"Y range for circle: ({lower_y:.2f}, {upper_y:.2f})")
    Y range for circle: (-1.00, 1.00)
    """
    from scipy.optimize import brentq, fsolve

    y_values = []
    x_samples = np.linspace(lower_x, upper_x, num_samples)
    min_y, max_y = y_search_range

    for x in x_samples:
        # Define function for this specific x
        def f_y(y):
            return implicit_func(x, y)

        # Try to find roots using brentq (requires bracketing interval)
        try:
            # First check if there's a sign change to bracket a root
            if f_y(min_y) * f_y(max_y) < 0:  # Sign change indicates a root
                root = brentq(f_y, min_y, max_y)
                y_values.append(root)
        except (ValueError, RuntimeError):
            pass

        # Also try fsolve for additional roots
        try:
            # Try multiple starting points to find different roots
            for y_start in np.linspace(min_y, max_y, 10):
                root = fsolve(f_y, y_start, full_output=True)
                if root[2] == 1:  # Solution found
                    y_val = root[0][0]
                    # Check if this is a valid solution and not already found
                    if abs(f_y(y_val)) < 1e-6 and not any(abs(y_val - yv) < 1e-3 for yv in y_values):
                        y_values.append(y_val)
        except:
            pass

    if not y_values:
        # If no solutions found, return the search range
        # return (min_y, max_y)
        return (-10, 10)

    y_values = np.array(y_values)
    lower_y = float(np.min(y_values))
    upper_y = float(np.max(y_values))

    return (lower_y, upper_y)
