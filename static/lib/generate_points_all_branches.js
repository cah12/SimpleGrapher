// generate_points_all_branches.js
// Port of Python `generate_points_all_branches` to JavaScript using mathjs.
// Exports `generatePointsAllBranches(expr, xMin, xMax, options)`.

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

  function defaultIfNull(v, d) {
    return v == null ? d : v;
  }

  function bisectRoot(f, a, b, tol = 1e-12, maxIter = 60) {
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

  function uniqueRounded(arr, digits = 12) {
    const seen = new Set();
    const out = [];
    for (const v of arr) {
      const key = Number(v).toFixed(digits);
      if (!seen.has(key)) {
        seen.add(key);
        out.push(Number(key));
      }
    }
    return out.sort((a, b) => a - b);
  }

  /**
   * Generate all real y(x) branches for implicit equation `expr(x,y)=0`.
   * @param {string|Object} expr - mathjs expression string or compiled node (must use variables `x` and `y`).
   * @param {number} xMin
   * @param {number} xMax
   * @param {Object} options
   * @returns {Array} branches - array of branches, each branch is [[x,y],...]
   */
  function generatePointsAllBranches(expr, xMin, xMax, options = {}) {
    const numX = defaultIfNull(options.numX, 400);
    let yMin = options.yMin;
    let yMax = options.yMax;
    const ySamples = defaultIfNull(options.ySamples, 400);
    const matchTol = options.matchTol; // may be null
    const fTol = defaultIfNull(options.fTol, 1e-15);
    const maxBisectionIter = defaultIfNull(options.maxBisectionIter, 60);
    const bisectionTol = defaultIfNull(options.bisectionTol, 1e-12);

    // Prepare expression
    let compiled;
    if (typeof expr === "string") {
      compiled = math.parse(expr).compile();
    } else if (expr && typeof expr.evaluate === "function") {
      compiled = expr;
    } else {
      throw new Error(
        "expr must be a mathjs expression string or compiled node",
      );
    }

    const f = function (xv, yv) {
      // attempt to evaluate, coerce to number
      try {
        const val = compiled.evaluate({ x: xv, y: yv });
        return Number(val);
      } catch (e) {
        return NaN;
      }
    };

    // Fallback y-range heuristic
    if (yMin == null || yMax == null) {
      const yGuess = Math.max(1.0, Math.abs(xMin), Math.abs(xMax)) * 10.0;
      if (yMin == null) yMin = -yGuess;
      if (yMax == null) yMax = yGuess;
    }

    const matchTolerance =
      matchTol == null ? Math.max(0.1, (yMax - yMin) * 0.1) : matchTol;

    const xVals = linspace(xMin, xMax, numX);
    const yGrid = linspace(yMin, yMax, ySamples);

    let branches = [];
    let prevRoots = [];

    for (let xi of xVals) {
      // evaluate f on yGrid
      const fvals = new Array(yGrid.length);
      let failed = false;
      for (let i = 0; i < yGrid.length; i++) {
        const v = f(xi, yGrid[i]);
        if (!isFinite(v) && !Number.isNaN(v)) {
          // Infinity
          fvals[i] = v;
        } else if (Number.isNaN(v)) {
          fvals[i] = NaN;
        } else {
          fvals[i] = v;
        }
      }

      const rootsThisX = [];

      // near-zero grid points
      for (let i = 0; i < fvals.length; i++) {
        const fv = fvals[i];
        if (!Number.isNaN(fv) && Math.abs(fv) <= fTol) {
          rootsThisX.push(yGrid[i]);
        }
      }

      // sign-change brackets
      for (let i = 0; i < fvals.length - 1; i++) {
        const aVal = fvals[i];
        const bVal = fvals[i + 1];
        if (Number.isNaN(aVal) || Number.isNaN(bVal)) continue;
        if (!isFinite(aVal) || !isFinite(bVal)) continue;
        if (aVal * bVal < 0) {
          const a = yGrid[i],
            b = yGrid[i + 1];
          try {
            const root = bisectRoot(
              (y) => f(xi, y),
              a,
              b,
              bisectionTol,
              maxBisectionIter,
            );
            if (root != null && !Number.isNaN(root)) rootsThisX.push(root);
          } catch (e) {
            // ignore
          }
        }
      }

      if (rootsThisX.length === 0) {
        prevRoots = [];
        continue;
      }

      const roots = uniqueRounded(rootsThisX, 12);
      const assigned = new Array(roots.length).fill(false);
      const newPrevRoots = [];

      // match to existing branches
      for (let bi = 0; bi < prevRoots.length; bi++) {
        const brow = prevRoots[bi];
        if (roots.length === 0) {
          newPrevRoots.push(brow);
          continue;
        }
        let diffs = roots.map((r) => Math.abs(r - brow));
        let bestIdx = diffs.indexOf(Math.min(...diffs));
        if (
          bestIdx >= 0 &&
          diffs[bestIdx] <= matchTolerance &&
          !assigned[bestIdx]
        ) {
          // append to branch
          branches[bi].push([Number(xi), Number(roots[bestIdx])]);
          assigned[bestIdx] = true;
          newPrevRoots.push(roots[bestIdx]);
        } else {
          // branch disappears at this x
          newPrevRoots.push(brow);
        }
      }

      // new branches for unassigned roots
      for (let ri = 0; ri < roots.length; ri++) {
        if (!assigned[ri]) {
          branches.push([[Number(xi), Number(roots[ri])]]);
          newPrevRoots.push(roots[ri]);
        }
      }

      prevRoots = newPrevRoots;
    }

    // keep only branches with at least 2 points
    branches = branches.filter((br) => br.length >= 2);

    return branches;
  }

  // Export
  if (typeof module !== "undefined" && module.exports) {
    module.exports = { generatePointsAllBranches };
  } else {
    global.generatePointsAllBranches = generatePointsAllBranches;
  }
})(this);

/*
Example usage (in browser with mathjs loaded or in Node with `mathjs`):

const expr = 'y^2 - x';
const branches = generatePointsAllBranches(expr, -1, 10, { numX: 200, ySamples: 400 });
console.log(branches.length, branches[0] && branches[0].slice(0,5));

*/
