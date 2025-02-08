// console.log("object");
async function solve_for(exp, c) {
  let res = [];

  //uncomment the code below to test in the in the small UI
  /*try {
    const _res = await axios.post("/csolve", { exp: exp, var: c });
    res = _res.data.result;
    if (!res.length) {
      $("#out").html(`Unable to find a solution`);
    } else {
      $("#out").html(`${c} = ${res.toString()}`);
    }
  } catch (error) {
    console.log(error.toString());
    $("#out").html(error.toString());
  }*/

  try {
    res = await axios.post("/csolve", { exp: exp, var: c });
  } catch (error) {
    console.error(error);
  }
  //console.log(res.data.result);
  return res.data.result;
}
