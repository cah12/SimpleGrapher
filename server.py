from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

import sympy as sp
from sympy.calculus.util import continuous_domain

# The Degree To Radian Code
#--------------------------#

mode_deg_rad = "deg"

# def setMode(mode):
#     mode_deg_rad = mode

def deg2rad(d):
    return sp.N(d*sp.pi/180)

original_sin = sp.sin
def sin(d):
    if mode_deg_rad == "deg":
        # print("deg")
        return original_sin(deg2rad(d))
    # print("rad")
    return original_sin(d)
sp.sin = sin

original_cos = sp.cos
def cos(d):
    if mode_deg_rad == "deg":
        # print("deg")
        return original_cos(deg2rad(d))
    # print("rad")
    return original_cos(d)
sp.cos = cos

original_tan = sp.tan
def tan(d):
    if mode_deg_rad == "deg":
        # print("deg")
        return original_tan(deg2rad(d))
    # print("rad")
    return original_tan(d)
sp.tan = tan

original_sec = sp.sec
def sec(d):
    if mode_deg_rad == "deg":
        # print("deg")
        return 1/original_cos(deg2rad(d))
    # print("rad")
    return 1/original_cos(d)
sp.sec = sec

original_csc = sp.csc
def csc(d):
    if mode_deg_rad == "deg":
        # print("deg")
        return 1/original_sin(deg2rad(d))
    # print("rad")
    return 1/original_sin(d)
sp.csc = csc

original_cot = sp.cot
def cot(d):
    if mode_deg_rad == "deg":
        # print("deg")
        return 1/original_tan(deg2rad(d))
    # print("rad")
    return 1/original_tan(d)
sp.cot = cot




def solve_for(exp, c):
    # print(exp)     
    while exp.find("^") != -1:
        exp = exp.replace("^", "**")
    
    v = sp.Symbol(c)
    arr = exp.split('=')

    
    
    if len(arr) == 1:
        arr.append("0")  
    
    result = []
    try:
        if arr[1] == "0":
            solutions = sp.solveset(sp.parse_expr(arr[0]), v)        
        elif arr[1] != "0":
            # print(arr[0])
            # print(arr[1])
            # print(v)
            s = sp.solveset(sp.Eq(sp.parse_expr(arr[0]), sp.parse_expr(arr[1])), v)
            # print(s)
            # print(type(s)==sp.ConditionSet)
            solutions = list(s) 
            # print(solutions)
                  
        for solution in solutions:
            num, denom = solution.as_numer_denom()
            if num.is_polynomial():
                num = sp.factor(num)  

            if denom.is_polynomial():
                denom = sp.factor(denom)             

            solution = (num)/(denom) 
            solution = solution.simplify()
            
            _str = str(solution)
            while _str.find("**") != -1:
                _str = _str.replace("**", "^")
            if _str.find("Piecewise") == -1 and _str.find("I") == -1:                
                result.append(_str)                               
        
    except BaseException as error:
        print(error)
        result = []    
    
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
    


def discontinuities(exp_, lower, upper, _var):  
    discount = []  
    v = sp.Symbol(_var)

    _exp=sp.parse_expr(exp_) 
    num, denom = _exp.as_numer_denom()    
    if num.is_polynomial():
        num = sp.factor(num)  

    if denom.is_polynomial():
        denom = sp.factor(denom)             

    solution = (num)/(denom) 
    solution = solution.factor()
    num, denom = solution.as_numer_denom()
    ds = sp.solveset(denom, v, sp.Interval(lower, upper, left_open=True, right_open=True))
    # print(ds)
    # print(type(ds))
    # print(type(ds)==sp.FiniteSet)
    

    if type(ds)==sp.FiniteSet:
        for sol in list(ds):
            discount.append(float(sol.evalf()))   

    # print(_exp)     

    if len(discount) == 0:
        # d = sp.factor(_exp) 
        d= list(sp.singularities(_exp, v))
        for sol in d:
            discount.append(float(sol.evalf()))

    discount = list(set(discount))
    discount.sort()   
    
    return discount
    
   

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
    # print("solve")
    return jsonify({"result": result})    


@app.route("/mode", methods=['POST'])
def mode():
    global mode_deg_rad
    data = request.get_json()  
    m = data["mode"]   
    mode_deg_rad = m
    return jsonify({"mode": m})   


@app.route("/discontinuity", methods=['POST'])
def discontinuity():
    data = request.get_json()    
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    while _exp.find("^") != -1:
        _exp = _exp.replace("^", "**")
    
    # print("discontinuity")
    discont = discontinuities(_exp, lower, upper, _var) 

    return jsonify({"discontinuities": discont})  


@app.route("/points", methods=['POST'])
def points():
    data = request.get_json()    
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    while _exp.find("^") != -1:
        _exp = _exp.replace("^", "**")
    
    discont = discontinuities(_exp, lower, upper, _var) 
    inflectn_points = inflection_points(_exp, lower, upper, _var) 
    turn_points = turning_points(_exp, lower, upper, _var)     
    return jsonify({
        "discontinuities": discont,
        "inflection_points": inflectn_points,
        "turning_points": turn_points
        })      
    

if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=3500, threads=100)



