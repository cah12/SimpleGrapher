
import numpy as np
import sympy as sp

def custom_sin_mode(arg):
    return np.sin(np.deg2rad(arg))


def custom_cos_mode(arg):
    return np.cos(np.deg2rad(arg))


def custom_tan_mode(arg):
    return np.tan(np.deg2rad(arg))


def custom_cot_mode(arg):
    return 1 / np.tan(np.deg2rad(arg))


def custom_sec_mode(arg):
    return 1 / np.cos(np.deg2rad(arg))


def custom_csc_mode(arg):
    return 1 / np.sin(np.deg2rad(arg))


def custom_asin_mode(arg):
    return np.rad2deg(np.arcsin(arg))


def custom_acos_mode(arg):
    return np.rad2deg(np.arccos(arg))


def custom_atan_mode(arg):
    return np.rad2deg(np.arctan(arg))


def custom_acot_mode(arg):
    return np.rad2deg(np.arccot(arg))


def custom_asec_mode(arg):
    return np.rad2deg(np.arccos(1/arg))


def custom_acsc_mode(arg):
    return np.rad2deg(np.arcsin(1/arg))


custom = {
    "sin_mode": custom_sin_mode,
    "cos_mode": custom_cos_mode,
    "tan_mode": custom_tan_mode,
    "cot_mode": custom_cot_mode,
    "sec_mode": custom_sec_mode,
    "csc_mode": custom_csc_mode,
    "asin_mode": custom_asin_mode,
    "acos_mode": custom_acos_mode,
    "atan_mode": custom_atan_mode,
    "acot_mode": custom_acot_mode,
    "asec_mode": custom_asec_mode,
    "acsc_mode": custom_acsc_mode
}