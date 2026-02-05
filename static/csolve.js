// console.log("object");
async function solve_for(exp, c) {
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
    res = await axios.post("/csolve", { exp: exp, var: c });
  } catch (error) {
    console.error(error);
  }
  // console.log(typeof res.data.result[0]);
  return res.data.result;
}

async function numeric(exp, lower, upper, indepVar, numOfPoints = 200) {
  try {
    res = await axios.post("/numeric", {
      exp: exp,
      lower: lower,
      upper: upper,
      var: indepVar,
      numOfPoints: numOfPoints,
    });
    // console.log(res.data);
    return res.data;
  } catch (error) {
    console.error(error);
    return [];
  }
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

async function discontinuity(exp, lower, upper, indepVar) {
  try {
    res = await axios.post("/discontinuity", {
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
