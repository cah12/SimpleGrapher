from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

import sympy as sp
import numpy as np

# The Degree To Radian Code
#--------------------------#
####################################################
###########################################################
mode_deg_rad = "deg"



def deg2rad(d):
    return sp.N(d*sp.pi/180)

original_sin = sp.sin
def sin(d):
    if mode_deg_rad == "deg":
        return original_sin(deg2rad(d))
    return original_sin(d)
sp.sin = sin

original_cos = sp.cos
def cos(d):
    if mode_deg_rad == "deg":
        return original_cos(deg2rad(d))
    return original_cos(d)
sp.cos = cos

original_tan = sp.tan
def tan(d):
    if mode_deg_rad == "deg":
        return original_tan(deg2rad(d))
    return original_tan(d)
sp.tan = tan

original_sec = sp.sec
def sec(d):
    if mode_deg_rad == "deg":
        return 1/original_cos(deg2rad(d))
    return 1/original_cos(d)
sp.sec = sec

original_csc = sp.csc
def csc(d):
    if mode_deg_rad == "deg":
        return 1/original_sin(deg2rad(d))
    return 1/original_sin(d)
sp.csc = csc

original_cot = sp.cot
def cot(d):
    if mode_deg_rad == "deg":
        return 1/original_tan(deg2rad(d))
    return 1/original_tan(d)
sp.cot = cot


def pyExpToJsExp(s):
    while s.find("**") != -1:
        s = s.replace("**", "^")
    return s


def jsExpToPyExp(s):
    while s.find("^") != -1:
        s = s.replace("^", "**")
    return s

def solve_for(exp, c):
    # while exp.find("^") != -1:
    #     exp = exp.replace("^", "**")

    exp = jsExpToPyExp(exp)
    
    v = sp.Symbol(c)
    arr = exp.split('=')

    if len(arr) == 2:
        if arr[0] == c:
            return[pyExpToJsExp(arr[1])]
        if arr[1] == c:
            return[pyExpToJsExp(arr[0])]

    
    if len(arr) == 1:
        arr.append("0")  
    
    result = []
    try:
        if arr[1] == "0":
            solutions = sp.solveset(sp.parse_expr(arr[0]), v)        
        else:
            s = sp.solveset(sp.Eq(sp.parse_expr(arr[0]), sp.parse_expr(arr[1])), v)
            solutions = list(s)            

        l = len(solutions)

        if l > 2:  
            if l % 2 == 0:
                solutions = [solutions[0], solutions[1]] 
            else:
                solutions = [solutions[0]] 

        for solution in solutions:
            num, denom = solution.as_numer_denom()
            if denom != 1:                
                if num.is_polynomial():
                    num = sp.factor(num)  

                if denom.is_polynomial():
                    denom = sp.factor(denom)             

                solution = (num)/(denom) 
                solution = solution.simplify()
            
            _str = str(solution)
            # while _str.find("**") != -1:
            #     _str = _str.replace("**", "^")
            _str = pyExpToJsExp(_str)
            if _str.find("Piecewise") == -1 and _str.find("I") == -1:                
                result.append(_str)                               
        
    except BaseException as error:
        print(error)
        return []     

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

    if denom == 1:
        return []  
    
    if num.is_polynomial():
        num = sp.factor(num)  

    if denom.is_polynomial():
        denom = sp.factor(denom)             

    
    solution = sp.factor(num/denom)
    num, denom = solution.as_numer_denom()
    ds = sp.solveset(denom, v, sp.Interval(lower, upper, left_open=True, right_open=True))   
    

    if type(ds)==sp.FiniteSet:
        for sol in list(ds):
            try:
                v = float(sol.evalf())
                discount.append(v)
            except:
                pass
  
    
    if len(discount) == 0:
        try:
            d= list(sp.singularities(_exp, v))
            for sol in d:
                try:
                    v = float(sol.evalf())
                    discount.append(v)
                except:
                    pass

            if len(discount):
                arr = []
                for s in discount:
                    if s < upper and s > lower:
                        arr.append(s)
                discount = arr
        except:
            pass

            

    
    discount.sort()   
    discount = list(map(float, discount))
    
    return discount   
 
#Checks if two sympy expressions are approximately equal by evaluating them numerically.
def approx_equal(expr1, expr2, tol=1e-9):    
    # Identify all free symbols in both expressions
    free_symbols = set(expr1.free_symbols) | set(expr2.free_symbols)
    
    # If there are symbols, substitute some random numerical values
    if free_symbols:
        # Generate random values for the symbols
        # A simple approach for demonstration; a more robust solution might use
        # a range of values or multiple test cases.
        values = {symbol: np.random.uniform(-10, 10) for symbol in free_symbols}
        
        # Evaluate expressions numerically
        num_expr1 = float(expr1.subs(values).evalf())
        num_expr2 = float(expr2.subs(values).evalf())
    else:
        # If no free symbols, just evaluate the expressions
        num_expr1 = float(expr1.evalf())
        num_expr2 = float(expr2.evalf())
        
    # Use numpy.allclose for a robust floating-point comparison
    return np.allclose(num_expr1, num_expr2, atol=tol)
   

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

    if _exp == None:
        return []

    # while _exp.find("^") != -1:
    #     _exp = _exp.replace("^", "**")

    _exp = jsExpToPyExp(_exp)
    
    # print("discontinuity")
    discont = discontinuities(_exp, lower, upper, _var) 
    # print(discont)

    return jsonify({"discontinuities": discont})  


@app.route("/points", methods=['POST'])
def points():
    data = request.get_json()    
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    # while _exp.find("^") != -1:
    #     _exp = _exp.replace("^", "**")

    _exp = jsExpToPyExp(_exp)
    
    discont = discontinuities(_exp, lower, upper, _var) 
    inflectn_points = inflection_points(_exp, lower, upper, _var) 
    turn_points = turning_points(_exp, lower, upper, _var)     
    return jsonify({
        "discontinuities": discont,
        "inflection_points": inflectn_points,
        "turning_points": turn_points
        })      
    

if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=8080, threads=100)



