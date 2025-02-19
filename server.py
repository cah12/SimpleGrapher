from flask import Flask, render_template, request, jsonify, make_response
from waitress import serve

import sympy as sp
from sympy.calculus.util import continuous_domain

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
    v = sp.Symbol(c)
    arr = exp.split('=')
    
    if len(arr) == 1:
        arr.append("0")  
    
    result = []
    try:
        if arr[1] == "0":
            solutions = sp.solveset(sp.parse_expr(arr[0]), v) 

        elif arr[1] != "0":
            solutions = list(sp.solveset(sp.Eq(sp.parse_expr(arr[0]), sp.parse_expr(arr[1])), v))    
        
        for solution in solutions:
            num, denom = solution.as_numer_denom()
            # print(num)
            # print(denom)
            if num.is_polynomial():
                # print(456)
                num = sp.factor(num)  

            if denom.is_polynomial():
                # print(556)
                denom = sp.factor(denom)             

            # print(num)
            # print(denom)
            solution = (num)/(denom) 
            # print(solution) 
            solution = solution.simplify()
            # print(solution)    
            
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

    # result = list(set(result)) 
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
    


def discontinuities(_exp, lower, upper, _var):
    try:
        right_bound = None
        discont = []
        domain = []
        fn = sp.parse_expr(_exp)
        u = continuous_domain(fn, sp.Symbol(_var), sp.Interval(lower, upper))
        
        
        if isinstance(u,sp.Interval):
            domain.append(list(map(float, u)))
            # return {"discont":discont, "domain":domain}
        
        if isinstance(u,sp.Union):
            for subset in u.args:            
                d = list(map(float, subset.boundary))
                domain.append(d)
                if right_bound != None and right_bound==d[0]:
                    discont.append(d[0])
                right_bound = d[1]            
            # return {"discont":discont, "domain":domain}         
        
    except:
        pass
    
    return {"discont":discont, "domain":domain} 

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
    _exp = data["exp"]
    _var = data["var"]
    lower = data["lower"]
    upper = data["upper"]

    while _exp.find("^") != -1:
        _exp = _exp.replace("^", "**")
    
    discont_domain = discontinuities(_exp, lower, upper, _var) 
    inflectn_points = inflection_points(_exp, lower, upper, _var) 
    turn_points = turning_points(_exp, lower, upper, _var)     
    return jsonify({
        "discontinuities": discont_domain["discont"],
        "domain": discont_domain["domain"],
        "inflection_points": inflectn_points,
        "turning_points": turn_points
        })      
    

if __name__ == "__main__":
    serve(app, host="0.0.0.0", port=3500, threads=100)



