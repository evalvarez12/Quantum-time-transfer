# -*- coding: utf-8 -*-


LOW_COST = 10.
MID_COST = 150.
HIGH_COST = 400.

def weight(a, b, c, d):
    return "calculated weight of structure"

def frequency(a, b, c, d):
    return "calculated resonant frequency"

def freq_penalty(freq):
    # Example linear piecewise penalty function -
    #   increasing cost for frequencies below 205 or above 395
    if freq < 205:
        return MID_COST * (205 - freq)
    elif freq < 395:
        return 0.
    else:
        return MID_COST * (freq - 395)

def stress_fraction(a, b, c, d):
    return "calculated stress / failure criteria"

def stress_penalty(stress_frac):
    # Example linear piecewise penalty function -
    #   low extra cost for stress fraction below 0.85,
    #   high extra cost for stress fraction over 0.98
    if stress_frac < 0.85:
        return LOW_COST * (0.85 - stress_frac)
    elif stress_frac < 0.98:
        return 0.
    else:
        return HIGH_COST * (stress_frac - 0.98)

def overall_fitness(parameter_vector):
    a, b, c, d = parameter_vector
    return (
        # D'oh! it took me a while to get this right -
        # we want _minimum_ weight and _minimum_ penalty
        # to get _maximum_ fitness.
       -weight(a, b, c, d)
      - freq_penalty(frequency(a, b, c, d))
      - stress_penalty(stress_fraction(a, b, c, d)
    )


from scipy.optimize import fmin

initial_guess = [29., 45., 8., 0.06]
result = fmin(lambda x: -overall_fitness(x), initial_guess, maxfun=100000, full_output=True)

