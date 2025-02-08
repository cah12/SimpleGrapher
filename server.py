from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

from sympy import *


def solve_for(exp, c):
    # print(exp)
    while exp.find("^") != -1:
        exp = exp.replace("^", "**")

    # print(exp)
    var = Symbol(c)
    arr = exp.split('=')
    
    if len(arr) == 1:
        arr.append("0")  
    
    result = []
    try:
        if arr[1] == "0":
            solutions = solve(parse_expr(arr[0]), var) 

        elif arr[1] != "0":
            solutions = solve(Eq(parse_expr(arr[0]), parse_expr(arr[1])), var)    
        # print(solutions)    
        for solution in solutions:
            _str = str(solution)
            while _str.find("**") != -1:
                _str = _str.replace("**", "^")
            if _str.find("Piecewise") == -1 and _str.find("I") == -1:                
                result.append(_str)
                #print(result)                
        
    except BaseException as error:
        #print('An exception occurred: {}'.format(error))
        # result = 'An exception occurred: {}'.format(error)
        result = []
    
    # print(result)

    result = list(set(result))
    # print("")
    # print("")
    # print("")
    # print(result)
    return result




app = Flask(__name__)
@app.route("/")
@app.route("/index")
def index():
    return render_template("index.html")

@app.route("/csolve", methods=['POST'])
def csolve():
    data = request.get_json()    
    exp = data["exp"]
    var = data["var"]
    # print(exp)
    # print(var)
    result = solve_for(exp, var) 
    # print(result)   
    return jsonify({"result": result})    


if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=3500)



