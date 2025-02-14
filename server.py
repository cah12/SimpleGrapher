from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

import sympy as sp

# The Degree To Radian Code
#--------------------------#

# deg_mode = "deg"

# def deg2rad(d):
#     return sp.N(d*sp.pi/180)

# _sin = sp.sin
# def sin(d):
#     if deg_mode == "deg":
#         return _sin(deg2rad(d))
#     return _sin(d)

# sp.sin = sin

# _cos = sp.cos
# def cos(d):
#     if deg_mode == "deg":
#         return _cos(deg2rad(d))
#     return _cos(d)

# sp.cos = cos




def solve_for(exp, c):
    # print(exp)
    while exp.find("^") != -1:
        exp = exp.replace("^", "**")

    # print(exp)
    var = sp.Symbol(c)
    arr = exp.split('=')
    
    if len(arr) == 1:
        arr.append("0")  
    
    result = []
    try:
        if arr[1] == "0":
            solutions = sp.solve(sp.parse_expr(arr[0]), var) 

        elif arr[1] != "0":
            solutions = sp.solve(sp.Eq(sp.parse_expr(arr[0]), sp.parse_expr(arr[1])), var)    
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


def inflection_points(expr, lower, upper, var):
    return [] 


def turning_points(expr, lower, upper, var):
    try:
        d = sp.diff(sp.parse_expr(expr), sp.Symbol(var))
        dis = sp.solveset(d, var, sp.Interval(lower, upper, left_open=True, right_open=True))        
        return list(map(float, dis))        
    except:
        return []  
    

def inflection_points(expr, lower, upper, var):
    try:
        x = sp.Symbol(var)
        d = sp.diff(sp.parse_expr(expr), x, x)
        dis = sp.solveset(d, var, sp.Interval(lower, upper))        
        return list(map(float, dis))        
    except:
        return []  



# def discontinuities(_exp, lower, upper, var):
#     try:
#         fn = sp.parse_expr(_exp)
#         n, d = fn.as_numer_denom()  
#         # print(d)
#         dis = sp.solveset(d, var, sp.Interval(lower, upper, left_open=True, right_open=True))
#         ls = list(map(float, dis))
#         # print(ls)
#         ls.sort() 
#         return  ls      
#     except:
#         return []

# def discontinuities(_exp, lower, upper, var):
#     try:
#         res = []
#         fn = sp.parse_expr(_exp) 
#         # print("hello")
#         ts = list(sp.Add.make_args(fn))
#         print(ts)
#         for t in ts:
#             n, dns = t.as_numer_denom() 
#             dis = sp.solve(dns, var)
            
#             for d in dis:
#                 _d = float(d)
#                 if _d < upper and _d > lower:
#                     try:
#                         res.index(_d)
#                     except:                    
#                         res.append(float(d))            
        
#         return  res      
#     except:
#         print("error")
#         return []
    
    

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
    result = solve_for(exp, var) 
    return jsonify({"result": result})    


# @app.route("/mode", methods=['POST'])
# def mode():
#     data = request.get_json()  
#     m = data["mode"]
#     mode_deg_rad = m
    
#     return jsonify({"mode": mode_deg_rad_})    


@app.route("/points", methods=['POST'])
def points():
    data = request.get_json()    
    exp = data["exp"]
    var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    while exp.find("^") != -1:
        exp = exp.replace("^", "**")
    
    discont = discontinuities(exp, lower, upper, var) 
    inflectn_points = inflection_points(exp, lower, upper, var) 
    turn_points = turning_points(exp, lower, upper, var)     
    return jsonify({
        "discontinuities": discont,
        "inflection_points": inflectn_points,
        "turning_points": turn_points
        })      
    

if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=3500)



