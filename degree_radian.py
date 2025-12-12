import sympy as sp

mode_deg_rad = "deg"


class sin_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        # Optional: define evaluation logic, e.g. for specific numerical inputs
        if mode_deg_rad == "deg":
            arg = arg * sp.pi / 180
        return sp.sin(arg)
        # otherwise return as a symbolic MySin
        pass


class cos_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        if mode_deg_rad == "deg":
            arg = arg * sp.pi / 180
        return sp.cos(arg)
        pass


class tan_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        if mode_deg_rad == "deg":
            arg = arg * sp.pi / 180
        return sp.tan(arg)
        pass


class cot_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        if mode_deg_rad == "deg":
            arg = arg * sp.pi / 180
        return sp.cot(arg)
        pass


class sec_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        if mode_deg_rad == "deg":
            arg = arg * sp.pi / 180
        return sp.sec(arg)
        pass


class csc_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        if mode_deg_rad == "deg":
            arg = arg * sp.pi / 180
        return sp.csc(arg)
        pass


class asin_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        v = sp.asin(arg)
        if mode_deg_rad == "deg":
            v = v * 180 / sp.pi
        return v
        pass


class acos_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        v = sp.acos(arg)
        if mode_deg_rad == "deg":
            v = v * 180 / sp.pi
        return v
        pass


class atan_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        v = sp.atan(arg)
        if mode_deg_rad == "deg":
            v = v * 180 / sp.pi
        return v
        pass


class acot_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        v = sp.acot(arg)
        if mode_deg_rad == "deg":
            v = v * 180 / sp.pi
        return v
        pass


class asec_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        v = sp.asec(arg)
        if mode_deg_rad == "deg":
            v = v * 180 / sp.pi
        return v
        pass


class acsc_mode(sp.Function):
    @classmethod
    def eval(cls, arg):
        v = sp.acsc(arg)
        if mode_deg_rad == "deg":
            v = v * 180 / sp.pi
        return v
        pass


def trig_substitutions(expr):
    if mode_deg_rad == "rad":
        return expr
    with sp.evaluate(False):
        expr = expr.subs(sp.sin, sin_mode)
        expr = expr.subs(sp.cos, cos_mode)
        expr = expr.subs(sp.tan, tan_mode)
        expr = expr.subs(sp.cot, cot_mode)
        expr = expr.subs(sp.sec, sec_mode)
        expr = expr.subs(sp.csc, csc_mode)
        expr = expr.subs(sp.asin, asin_mode)
        expr = expr.subs(sp.acos, acos_mode)
        expr = expr.subs(sp.atan, atan_mode)
        expr = expr.subs(sp.acot, acot_mode)
        expr = expr.subs(sp.asec, asec_mode)
        expr = expr.subs(sp.acsc, acsc_mode)

    return expr


def set_mode(new_mode):
    global mode_deg_rad
    if new_mode in ["deg", "rad"]:
        mode_deg_rad = new_mode
