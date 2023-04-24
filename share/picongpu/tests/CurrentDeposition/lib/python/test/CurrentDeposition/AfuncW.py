import numpy as np


# functions for arrays (_a)
def NGP_a(x):
    """generates hat function with width 1 (symmetrical around x=0) -
    assignement function of Order 0
    """
    y = np.empty(len(x))

    for i in range(len(x)):
        if abs(x[i]) < 1/2:
            y[i] = 1
        elif abs(x[i]) == 1/2:
            y[i] = 1/2
        else:
            y[i] = 0
    return y


def CIC_a(x):
    """generates hat triangle with width 2 (symmetrical around x=0) -
    assignement function of Order 1
    """
    y = np.empty(len(x))

    for i in range(len(x)):
        if abs(x[i]) < 1:
            y[i] = 1 - abs(x[i])
        else:
            y[i] = 0
    return y


def TSC_a(x):
    """assignementfunction of order 2"""
    y = np.empty(len(x))

    for i in range(len(x)):
        if abs(x[i]) < 3/2 and abs(x[i]) >= 1/2:
            y[i] = ((3/2 - abs(x[i]))**2)/2
        elif abs(x[i]) < 1/2:
            y[i] = (-x[i]**2 + 3/4)
        else:
            y[i] = 0
    return y


# functions for single values
def NGP(x):
    """generates hat function with width 1 (symmetrical around x=0) -
    assignement function of Order 0
    """

    if abs(x) < 1/2:
        y = 1
    elif abs(x) == 1/2:
        y = 1/2
    else:
        y = 0
    return y


def CIC(x):
    """generates hat triangle with width 2 (symmetrical around x=0) -
    assignement function of order 1 for single values
    """

    if abs(x) < 1:
        y = 1 - abs(x)
    else:
        y = 0
    return y


def TSC(x):
    """assignementfunction of order 2" for single values"""

    if abs(x) < 3/2 and abs(x) >= 1/2:
        y = ((3/2 - abs(x))**2)/2
    elif abs(x) < 1/2:
        y = (-(x**2) + 3/4)
    else:
        y = 0
    return y


def PSQ(x):
    """assignementfunction of order 3 for single values"""

    if abs(x) < 2 and abs(x) >= 1:
        y = 1/6*(2 - abs(x))**3
    elif abs(x) < 1:
        y = (2/3 - x**2 + (abs(x)**3)/2)
    else:
        y = 0
    return y


def W(s1, s2, s3, s4, s5, s6):
    """Calculation of the current deposiotion; si are the values of the
    assignement function at the old and the new coordinates
    si = assignment-function(i - cell_position)"""

    W = 1/3*(s4*s5*s6 - s1*s5*s6 + s4*s2*s3 - s1*s2*s3) + (
        1/6*(s4*s2*s6 - s1*s2*s6 + s4*s5*s3 - s1*s5*s3))
    return W
