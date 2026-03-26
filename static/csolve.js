// console.log("object");
async function solve_for(math_mode, exp, c) {
  let res = [];

  //uncomment the code below to test in the in the small UI
  /* try {
    const _res = await axios.post("/csolve", { exp: exp, var: c });
    res = _res.data;
    if (!res) {
      $("#out").html(`Unable to find a solution`);
    } else {
      console.log(res);
      $("#out").html(JSON.stringify(res));
    }
  } catch (error) {
    console.log(error.toString());
    $("#out").html(error.toString());
  } */

  try {
    res = await axios.post("/csolve", {
      math_mode: math_mode,
      exp: exp,
      var: c,
    });
  } catch (error) {
    console.error(error);
  }
  // console.log(typeof res.data.result[0]);
  return res.data.result;
}

// Function to convert a Base64 string to an ArrayBuffer
function base64ToArrayBuffer(base64) {
  // Decode base64 to binary string
  const binaryString = atob(base64);
  const len = binaryString.length;
  const bytes = new Uint8Array(len);
  for (let i = 0; i < len; i++) {
    bytes[i] = binaryString.charCodeAt(i);
  }
  // Assume the numpy array is float32 and 2D, shape (N, 2)
  // Each pair is [x, y] as float32
  const floatArray = new Float64Array(bytes.buffer);
  const result = [];
  for (let i = 0; i < floatArray.length; i += 2) {
    // result.push([floatArray[i], floatArray[i + 1]]);
    let y = floatArray[i + 1];
    if (Math.abs(y) == 3.4e38) {
      y = 1e300 * Math.sign(y);
    }
    result.push(new Misc.Point(floatArray[i], y));
  }
  return result;
}

async function numeric(
  math_mode,
  exp,
  lower,
  upper,
  lowerY,
  upperY,
  autoScale,
  indepVar,
  numOfPoints = 200,
  has_discontinuity = false,
) {
  try {
    res = await axios.post(
      "/numeric",
      {
        exp: exp,
        lower: lower,
        upper: upper,
        lowerY: lowerY,
        upperY: upperY,
        autoScale: autoScale,
        var: indepVar,
        numOfPoints: numOfPoints,
        has_discontinuity: has_discontinuity,
        math_mode: math_mode,
      } /* ,
      { responseType: "arraybuffer" }, */,
    );
    const data = res.data;
    if (data.numpy_dtype !== "float32") {
      throw new Error("Unhandled numpy dtype:", data.numpy_dtype);
    }

    const array_of_branches = data.branches;
    const _array_of_branches = [];
    for (let i = 0; i < array_of_branches.length; i++) {
      _array_of_branches.push(base64ToArrayBuffer(array_of_branches[i]));
    }

    return {
      branches: _array_of_branches,
      discontinuities: data.discontinuities,
      large_range_span: data.large_range_span,
    };
  } catch (error) {
    console.error(error);
    return [];
  }
  // try {
  //   res = await axios.post("/numeric", {
  //     exp: exp,
  //     lower: lower,
  //     upper: upper,
  //     var: indepVar,
  //     numOfPoints: numOfPoints,
  //     has_discontinuity: has_discontinuity,
  //   });
  //   // console.log(res.data);
  //   return res.data;
  // } catch (error) {
  //   console.error(error);
  //   return [];
  // }
}

async function points(exp, lower, upper, indepVar) {
  try {
    res = await axios.post("/points", {
      exp: exp,
      lower: lower,
      upper: upper,
      var: indepVar,
    });
    // console.log(res.data);
    return res.data;
  } catch (error) {
    console.error(error);
    return [];
  }
}

async function turningPoints(exp, lower, upper, indepVar) {
  try {
    res = await axios.post("/turningPoints", {
      exp: exp,
      lower: lower,
      upper: upper,
      var: indepVar,
    });
    // console.log(res.data);
    return res.data;
  } catch (error) {
    console.error(error);
    return [];
  }
}

async function discontinuity(mode, exp, lower, upper, indepVar) {
  try {
    res = await axios.post("/discontinuity", {
      math_mode: mode,
      exp: exp,
      lower: lower,
      upper: upper,
      var: indepVar,
    });
    // console.log(res.data);
    return res.data;
  } catch (error) {
    console.error(error);
    return [];
  }
}

async function mode(_mode) {
  try {
    res = await axios.post("/mode", {
      mode: _mode,
    });
    // console.log(res.data);
    return res.data;
  } catch (error) {
    console.error(error);
    return null;
  }
}
