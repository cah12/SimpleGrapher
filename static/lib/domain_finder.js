// domain_finder.js
// Port of Python `closer_boundary` and `is_x_in_domain_numerical` to JavaScript using mathjs.
// Exports `closerBoundary(fn, currentBoundary, nextPoint, forward)`.

(function (global) {
  "use strict";

  // Ensure mathjs is available as `math` in the environment.
  const math = typeof require === "function" ? require("mathjs") : global.math;
  if (!math) {
    throw new Error(
      'mathjs is required (provide `math` global or use CommonJS require("mathjs")).',
    );
  }

  function linspace(a, b, n) {
    if (n === 1) return [a];
    const out = new Array(n);
    const step = (b - a) / (n - 1);
    for (let i = 0; i < n; i++) out[i] = a + step * i;
    return out;
  }

  function bisectRoot(f, a, b, tol = 1e-15, maxIter = 60) {
    let fa = f(a);
    let fb = f(b);
    if (!isFinite(fa) || !isFinite(fb)) return null;
    if (fa === 0) return a;
    if (fb === 0) return b;
    if (fa * fb > 0) return null; // no sign change

    let left = a,
      right = b,
      mid,
      fmid;
    for (let i = 0; i < maxIter; i++) {
      mid = 0.5 * (left + right);
      fmid = f(mid);
      if (!isFinite(fmid)) return null;
      if (Math.abs(fmid) <= tol || (right - left) / 2 < tol) return mid;
      if (fa * fmid <= 0) {
        right = mid;
        fb = fmid;
      } else {
        left = mid;
        fa = fmid;
      }
    }
    return 0.5 * (left + right);
  }

  /**
   * Check if a given x_val is in the domain by trying to find a valid real y.
   * @param {Function} fnXY - Function f(x, y) to check; should return a number.
   * @param {number} xVal - The x value to check.
   * @param {Array} yRange - [yMin, yMax] range to search for a valid y.
   * @param {number} tol - Tolerance for the root finder.
   * @returns {boolean} True if a valid y is found for the given x.
   */
  function isXInDomainNumerical(fnXY, xVal, yRange = [-100, 100], tol = 1e-15) {
    const scope = new Map();
    const fY = (y) => {
      try {
        const result = fnXY(xVal, y);
        // scope.set("x", xVal);
        // scope.set("y", y);
        // const result = math.evaluate(fnXY, scope);

        return Number(result);
      } catch (e) {
        return NaN;
      }
    };

    try {
      // Sample y range to find sign changes
      const ySamples = linspace(yRange[0], yRange[1], 100);
      const signs = ySamples.map((y) => Math.sign(fY(y)));

      // Look for sign changes
      for (let i = 0; i < signs.length - 1; i++) {
        if (signs[i] * signs[i + 1] < 0) {
          // Sign change found
          const result = bisectRoot(fY, ySamples[i], ySamples[i + 1], tol);
          if (result != null) {
            return true; // A valid y was found
          }
        }
      }

      // Check for near-zero values
      for (const y of ySamples) {
        if (Math.abs(fY(y)) < tol) {
          return true;
        }
      }

      return false;
    } catch (e) {
      return false;
    }
  }

  /**
   * Find a closer boundary point for the curve near current_boundary moving towards/away from next_point.
   * @param {Function} fn - Function f(x, y) returning a number.
   * @param {Array} currentBoundary - [x, y] of current boundary point.
   * @param {Array} nextPoint - [x, y] of next point.
   * @param {boolean} forward - If true, move forward; if false, move backward.
   * @returns {Array|null} [x, y] of closer boundary point or null if not found.
   */
  function closerBoundary(fn, currentBoundary, nextPoint, forward = true) {
    const node = math.parse(fn);
    const _fn = function (x, y) {
      scope.set("x", x);
      scope.set("y", y);
      return node.evaluate(scope);
    };
    const [cbx, cby] = currentBoundary;
    const [npx, npy] = nextPoint;

    const yStep = Math.abs(cby - npy);
    const xStep = Math.abs(cbx - npx);

    const stepFactor = 4;

    let yRange, xRange;
    if (forward) {
      yRange = [cby - stepFactor * yStep, cby + stepFactor * yStep];
      xRange = [cbx - stepFactor * xStep, cbx + stepFactor * xStep];
    } else {
      yRange = [cby + stepFactor * yStep, cby - stepFactor * yStep];
      xRange = [cbx + stepFactor * xStep, cbx - stepFactor * xStep];
      // Normalize ranges if backwards
      // yRange = [Math.min(yRange[0], yRange[1]), Math.max(yRange[0], yRange[1])];
      // xRange = [Math.min(xRange[0], xRange[1]), Math.max(xRange[0], xRange[1])];
    }

    const xSamples = linspace(xRange[0], xRange[1], 500);
    const ySamples = linspace(yRange[0], yRange[1], 1000);

    const scope = new Map();

    for (const x of xSamples) {
      if (isXInDomainNumerical(_fn, x, yRange)) {
        // Found an x in domain; now find the corresponding y
        for (let i = 0; i < ySamples.length - 1; i++) {
          const y0 = ySamples[i];
          const y1 = ySamples[i + 1];
          const f0 = _fn(x, y0);
          const f1 = _fn(x, y1);

          // Check for sign change
          if (!isNaN(f0) && !isNaN(f1) && isFinite(f0) && isFinite(f1)) {
            if (f0 * f1 < 0) {
              // Sign change found; use bisection to refine
              const root = bisectRoot((y) => _fn(x, y), y0, y1, 1e-12);
              if (root != null) {
                return [x, root];
              }
            }
          }
        }
      }
    }

    return null;
  }

  // Export
  if (typeof module !== "undefined" && module.exports) {
    module.exports = { closerBoundary, isXInDomainNumerical };
  } else {
    global.closerBoundary = closerBoundary;
    global.isXInDomainNumerical = isXInDomainNumerical;
  }
})(this);

/*
Example usage (in browser with mathjs loaded or in Node with `mathjs`):

// Define a function f(x, y) = y^6 + 8*y - x
const fn = (x, y) => Math.pow(y, 6) + 8*y - x;

const currentBoundary = [-6.784, -0.9276];
const nextPoint = [-6.683, -0.9034];

const result = closerBoundary(fn, currentBoundary, nextPoint, true);
console.log('Closer boundary point:', result);

*/
